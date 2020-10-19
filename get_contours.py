#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
    Copyright (C) 2020  Maria Oliver-Parera 
    <maria.oliver-parera@gipsa-lab.grenoble-inp.fr>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as published
    by the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.
"""

import sys
import cv2
import os
import math
import numpy as np
import argparse
from scipy.ndimage import gaussian_filter
from skimage.morphology import disk, opening

from flowio import readFlowFile

def normalizeImage(image, oldMin, oldMax, newMin, newMax):
    newIm = (image - oldMin)*(newMax - newMin)/(oldMax - oldMin) + newMin
    return newIm

def optical_strain(of, norm = 2, sigma = 1.8):
    u = of[...,0]
    v = of[...,1]
    
    if(sigma > 0):
        u = gaussian_filter(u, sigma)
        v = gaussian_filter(v, sigma)

    ex, uy = np.gradient(u)
    vx, ey = np.gradient(v)

    if norm == 1:
        e_xy = 0.5*(uy + vx)
        e_m2 = np.abs(ex) + 2*np.abs(e_xy) + np.abs(ey)
        e_m = np.sqrt(e_m2)
    elif norm == 2:
        e_xy = 0.5*(uy + vx)
        e_m2 = ex*ex + 2*e_xy*e_xy + ey*ey;
        e_m = np.sqrt(e_m2)
    elif norm == 3:
        e_m = np.abs(ex + ey)
    elif norm == 4:
        e_xy = 0.5*(uy + vx)
        e_m = np.abs(2*e_xy + ex*ey)
    elif norm == 5:
        e_xy = 0.5*(uy + vx)
        e_m = np.abs(ex + ey + 2*e_xy)
    else:
        print("norm not known")
        return -1
    return e_m

def get_otsu_threshold():
    return

def get_thresholds(T_influence, of_files, rho_l, rho_u, df, f_number, 
                   last_frame, nt, threshold_number, double_threshold, 
                   norm, sigma):
    
    of = readFlowFile(of_files[f_number])
    v_link = np.zeros_like(of)
    del of

    last_frame += df
    if last_frame >= nt:
        last_frame = nt
    
    epsilon_max = []
    epsilon_min = []
    while f_number < last_frame:
        of = readFlowFile(of_files[f_number])
        v_link = v_link + of
        e_m = optical_strain(v_link, norm, sigma);
        epsilon_max.append(np.max(e_m));
        epsilon_min.append(np.min(e_m));
        T_influence[f_number].append(threshold_number);
        f_number = f_number + 1;
        del of
    del e_m;
    maxmax = np.max(epsilon_max)
    minmin = np.min(epsilon_min)
    Tl = minmin + rho_l*(maxmax - minmin)
    Tu = maxmax;
    if double_threshold :
        Tu = maxmax - rho_u * (maxmax - minmin)
    
    del v_link
    return Tl, Tu, last_frame

def get_gaussian_weights(size, sigma):
    gaussian_weights = []
    middle = math.floor(size/2)
    for i in range(-middle, middle+1):
        if(size%2 == 0 and i == 0):
            continue
        weight = math.exp(-(i * i) / (2 * sigma * sigma));
        gaussian_weights.append(weight)

    return gaussian_weights

def weighted_threshold(thresholds, threshold_list, weight_type, shift):
    weighted_threshold = 0.0
    normalization = 0.0

    if(weight_type == 1):
        gaussian_weights = np.array(threshold_list.size())
        sigma = 1.0
        gaussian_weights = get_gaussian_weights(threshold_list.size(), sigma);
        count = 0;
        for thr in threshold_list: #"gaussian"
            weighted_threshold += gaussian_weights[count]*thr
            normalization += gaussian_weights[count]
            count = count + 1
        del gaussian_weights;
    elif(weight_type == 2):   #increasing
        count = 1
        for thr in thresholds:
            weighted_threshold += count*thr
            normalization += count
            count += shift
    elif(weight_type == 3): #decreasing
        count = threshold_list.size()
        for thr in thresholds:
            weighted_threshold = count*thr
            normalization += count
            count -= shift
    elif(weight_type == 4):   #no weights (monotonous)
        for thr in thresholds:
            weighted_threshold += thr
            normalization = normalization + 1

    return weighted_threshold/normalization

def apply_morphology(image):
    r = 2
    selem = disk(r)
    output = opening(image, selem)
    return output

def apply_threshold_frame(Tl, Tu, T_influence, of_files, df, f_number, 
                           last_frame, nt, weight_type, saveDir, norm, 
                           sigma, shift, threshold = True, morphology = True):
    real_tl = 0.0;

    of = readFlowFile(of_files[f_number])
    e_m = optical_strain(of, norm, sigma)
    output = e_m*255
    
    if threshold:
        e_m_threshold = np.zeros_like(e_m)
    
        if(len(T_influence[f_number]) == 1):
            real_tl = Tl[T_influence[f_number][0]]
        else:
            real_tl = weighted_threshold(Tl, T_influence[f_number], weight_type, shift)
        e_m_threshold[e_m > real_tl] = 255
        output = e_m_threshold
        if morphology:
            e_m_closing = apply_morphology(e_m_threshold)
            output = e_m_closing
    out_name = os.path.join(saveDir, of_files[f_number].split('/')[-1].split('.')[0] + '.png')
    cv2.imwrite(out_name, output)
    return

def get_contours(of_files, m, shift, rho_l, rho_u, weight_type, save_folder,
                 thr_type, norm, sigma, threshold = True, morphology = True):
    f_it = 0
    f_number = 0
    threshold_number = 0
    Tup = []
    Tlow = []
    nt = len(of_files)
    T_influence = [[] for _ in range(nt)]

    '''first sequence (we need it separately in case we are using the sliding window)'''
    double_threshold = False
    if thr_type == 2:
        double_threshold = True
    Tl, Tu, f_it = get_thresholds(T_influence, of_files, rho_l, rho_u, m, 
                                  f_number, f_it, nt, threshold_number, 
                                  double_threshold, norm, sigma)
    threshold_number = threshold_number + 1;
    Tlow.append(Tl)
    Tup.append(Tu)
    f_number += shift
    while (f_it < nt or f_number < nt):
        Tl, Tu, f_it = get_thresholds(T_influence, of_files, rho_l, rho_u, shift, 
                                      f_number, f_it, nt, threshold_number, 
                                      double_threshold, norm, sigma)
        threshold_number = threshold_number + 1
        Tlow.append(Tl)
        Tup.append(Tu)
        f_number += shift;

    ''' apply threshold frame by frame '''
    for i in range(0, nt):
        apply_threshold_frame(Tlow, Tup, T_influence, of_files, shift, i, f_it,
                              nt, weight_type, save_folder, norm, sigma, shift,
                              threshold, morphology)

    return

    

if __name__ == '__main__':
    print("Contour Detection of Moving Objects using Optical Strain")
    
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--files_dir", required=True, 
                        help="folder with the optical flow files")
    parser.add_argument("-o", "--save_dir", required=True, 
                        help="folder to save the optical strain files")
    parser.add_argument("-m", "--m", required=False, default = 10,
                        help="number of frames per window")
    parser.add_argument("-j", "--j", required=False, default = 5,
                        help="shift step: how many frames are jumped between windows")
    parser.add_argument("-s", "--sigma", required=False, default = 1.8,
                        help="mean of the gaussian applied at each optical flow")
    parser.add_argument("-n", "--norm", required=False, default = 2,
                        help="norm used for the opical strain")
    parser.add_argument("-t", "--thr_type", required=False, default = 1,
                        help="type of threshold applied")
    parser.add_argument("-l", "--rho_l", required=False, default = 0.05,
                        help="low percentile for the threshold of each frame")
    parser.add_argument("-u", "--rho_u", required=False, default = 0.0,
                        help="up percentile for the threshold of each frame")
    parser.add_argument("-w", "--weight_type", required=False, default = 2,
                        help="type of weight used to ponderated among all the" + 
                        "thresholds that affect one frame")
    parser.add_argument("-c", "--threshold", required=False, default = True,
                        help="boolean: True to apply threshold; False otherwise")
    parser.add_argument("-r", "--morphology", required=False, default = True,
                        help="boolean: True to apply morphology; False otherwise")
    args = vars(parser.parse_args())
    
    if not os.path.exists(args["save_dir"]):
        os.makedirs(args["save_dir"])
        
    of_names = sorted(os.listdir(args["files_dir"]))
    complete_names = [os.path.join(args["files_dir"], s) for s in of_names]    
    get_contours(complete_names, args["m"], args["j"], args["rho_l"], 
                 args["rho_u"], args["weight_type"], args["save_dir"], 
                 args["thr_type"], args["norm"], args["sigma"], 
                 args["threshold"], args["morphology"])
