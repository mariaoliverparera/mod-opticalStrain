-------------------------------
CONTOUR DETECTION OF MOVING OBJECTS USING OPTICAL STRAIN
-------------------------------

**Usage**

    get_contours.py [-h] -i FILES_DIR -o SAVE_DIR [-m M] [-j J] [-s SIGMA]
                       [-n NORM] [-t THR_TYPE] [-l RHO_L] [-u RHO_U]
                       [-w WEIGHT_TYPE] [-c THRESHOLD] [-r MORPHOLOGY]


**Required parameters**

    FILES_DIR : folder with the optical flow files
    SAVE_DIR : folder to save the optical strain files

**Optional parameters**

    M : number of frame per window (default : 10)
    J : number of frames jumped between windows (default : 5) corresponds to the euclidean norm)
    SIGMA : mean of the gaussian applied at each optical flow before applying the optical strain (default : 1.8)
    NORM : norm used to compute the optical strain (default : 2, which
    THR_TYPE : choice between only a lower threshold (thr_type = 1), which removed the low values; or an upper and lower threshold (thr_type = 2), which removes both, lower and upper values(default : 1)
    RHO_L : low percentile for the threshold of each frame, it is used in the computation of the thresholds (default : 0.05)
    RHO_U : up percentile for the threshold of each frame, it is used in the computation of the thresholds (default : 0.0)
    WEIGHT_TYPE : determines the weight associated to each threshold associated to each frame, After computing the threshold associated to each window, it is applied to each frame. But each frame can belong to more than one window, with the consequence of having more than one frame. (default : 2, which corresponds to an "increasing" ponderation, i.e, the first frame has less weight than the last one)
    THRESHOLD : return the values after a threshold of keep the values from the optical strain (default : True)
    MORPHOLOGY : apply a closing operator after threshold (default: True)

**Example**

    ./get_contours.py -i bear/ -o results

The images provided as example are obtained using the "faldoi" algorithm to compute optical flow (http://www.ipol.im/pub/art/2019/238/) over the sequence "bear" from DAVIS 2017 Unsupervised dataset (https://davischallenge.org/).

**Citation**
    If you use "get_contours", please cite the following paper:

    @article{olivercontour,
      title={Contour Detection of Multiple Moving Objects in Unconstrained Scenes using Optical Strain},
      author={Oliver-Parera, Maria and Muzeau, Julien and Ladret, Patricia and Bertolino, Pascal}
    }


**Copyright and license**

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as
    published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.

    You should have received a copy of the GNU Affero General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.

    -------
    LICENSE
    -------
    File flowio.py copied verbatim from its original sources.
