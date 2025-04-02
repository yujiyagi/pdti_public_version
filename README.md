# Potency Density Tensor Inversion (PDTI)

The advantages of introducing uncertainty into the Green's function in source inversion are described in Yagi & Fukahata (2011).<br>
The advantages of using a five-component basis double-couple are described in Shimizu, et al. (2020).<br>
The advantages of time-adaptive smoothing are described in Yamahsita et al. (2022).

## Build
```bash
$ git clone git@github.com:yujiyagi/pdti_public_version.git
$ cd PDTI
$ mkdir $HOME/data/   ! If there is no $HOME/data directory.
$ cp -a data/ak135_taup.*   $HOME/data/    ! copy P-wave tables
# We have prepared a P-wave tables calculated using TauP (https://www.seis.sc.edu/taup/), but we expect that you will create the tables yourself. 
$ cd src
$ make
$ cd 
$ echo "export PATH=$HOME/PathToPDTI/bin:PATH" >> $HOME/.bashrc
  ```

### Preparation for PDTI.
 "wave.obs" directory containing observed waveforms (Displacement waveform with corrected seismometer response)
 ```
 Example of waveform file name:  MAJOBHZ  (station code)+(comp)  
 Format of waveform file: 
   Tp(timing of P arrival, please set 10.0 sec), dt (sampling interval), nt (number of data point), lat, lon, tmax (timing of PP 
arrival), ndmax (Point of PP arrival), sig (Variance of observation error)
      x(i), i=1,nt
Example: 
    10.00   0.0250    15600    36.55   138.20     100.20      4008      0.330583
 -0.160335E-03
 -0.358657E-03
 -0.634570E-03
 -0.986220E-03
 -0.141550E-02
 ```
 structure.dat
  ```
Tp*  Ts*  n_h   ( vp_h(n)  vs_h(n)  density_h(n)  thinkness_h(n)  dumy  dumy, i=1,n_h)
n_s ( vp_s(n)  vs_s(n)  density_s(n)  thinkness_s(n), i=1,n_s)
n_i ( vp_s(n)  vs_s(n)  density_s(n)  thinkness_s(n), i=1,n_i)
Example:  Model AK135 (continental structure)
   1. 4. 3  5.8000  3.4600  2.4490  20.0  0.  0.  
            6.5000  3.8500  2.7142  15.0  0.  0. 
            8.0400  4.4800  3.2976   0.0  0.  0. 
         3  5.8000  3.4600  2.4490  20.0   
            6.5000  3.8500  2.7142  15.0   
            8.0400  4.4800  3.2976   0.0   
         1  5.8000  3.4600  2.4490   20.
 ```
epocenter.dat
  ```
    lat lon depth
Example:   22.013  95.922  15.0
 ```
station.list
  ```
     $stcd  $comp  $lat  $lon  $azimuth  $delta
Example:
TLY    BHZ   51.681   103.644     9.5    30.2  
TIXI   BHZ   71.634   128.867    12.4    53.1  
  ```
i_greenf
```
$nt  $dtg  $str  $dip  $rake $depth $vr 
$XX  $YY  $MN   $NN   $M0   $N0    $ICMN  $Rslip $nsurface

nt: number of points for green’s function
dtg:  sampling interval
str:   strike of model plane
dip:   dip of model plane
rake:   dummy
depth:   depth of hypocenter
vr:    assumed maximum rupture velocity
XX:   knot interval for strike direction
YY:   knot interval for dip direction
MN:  Total knot number for strike direction 
NN:  Total knot number for dip direction 
M0,N0:   Knot location for hypocenter 
ICMN:  Number of base double couples
Rslip:  dummy
nsurface: Setting the smoothing constraint conditions for the top edge (0, 1 or 2)
```
i_Cmatrix
```
$title 
$TR  $JTN  $vr $sift $p1 $p2 $p3 $st_max $r_s_t " 1 " $criterion 
$alpha  1.  1.   $nabic  $nsurface " 0  F "  1.  $gamma  $dump 
$stcd  $comp  $cp   $wight  $tl   $dt
TR: knot interval for potency-rate time function
JTN: maximum knot number for potency-rate time function
Sift, p1, p2, p3: dummy 
st_max: Total length of potency-rate time function 
r_s_t: resolution ratio in equation (24) of Yagi and Fukahata (2011)
criterion: Variable c in equation (10) of Yamaguchi et al. (2022)
alpha: hyperparameter alpha in equation (27) of Yagi and Fukahata (2011)
gamma: hyperparameter gamma in equation (27) of Yagi and Fukahata (2011)
nabic: dummy (set 1)
dump: Variable w in equation (32) of Yagi and Fukahata (2011)
cp: dummy (set 0)
wight: dummy (set 0.1)
tl: length of waveform 
dt: sampling interval for inversion
```

### Run PDTI
```bash
$ PDTI_get_knot_info < i_greenf     ! Calculating the position of spatial knots
$ PDTI_GreenPointSources < i_greenf  ! Calculation of GF
$ PDTI_get_ndj_main  < i_Cmatrix    ! Calculating the appropriate time window for each observation
$ PDTI_get_init_model_para_TmpSm < i_Cmatrix  !  Setting the initial value
$ PDTI_pre_inv_new  < i_Cmatrix
Repeat the following two lines of command until the change in the solution is sufficiently small. Please enter 1.0 for the alpha 
and gamma values.
$ PDTI_get_covariance_grn  < i_Cmatrix  ! Calculation of the covariance matrix of the GF error
$ PDTI_ABIC_TA < i_Cmatrix
When convergence is reached, check the values of alpha and gamma that minimise ABIC, and then modify i_Cmatrix.
$ PDTI_inversionA_TAd < i_Cmatrix
The results are output to fort.40. For details of this file, please refer to ‘sub.w_sol.f90’.
```

### reference
Please cite the following three papers.
  - [Yagi & Fukahata, 2011, GJI](https://doi.org/10.1111/j.1365-246X.2011.05043.x)  
  - [Shimizu, Yagi, Okuwaki & Fukahata, 2019, GJI](https://doi.org/10.1093/gji/ggz496) 
  - [Yamashita, Yagi, Okuwaki, Shimizu, Agata &  Fukahata, 2022, GJI](https://doi.org/10.1093/gji/ggac181)

If you calculate the Green function using GreenPointSources.f90, please cite [Kikuchi & Kanamori, 1991, BSSA](https://doi.org/10.1785/BSSA0810062335).
If you use the P-wave table calculated using TauP, please cite [Crotwell, Owens, Ritsema, 1999, SRL](https://doi.org/10.1785/gssrl.70.2.154).

