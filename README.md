# pSAR - Automated InSAR Processing Toolkit for Sentinel-1 TOPS Data

# Overview
pSAR is an automated InSAR (Interferometric Synthetic Aperture Radar) processing toolkit developed based on Python and GMTSAR. It is specifically designed for Sentinel-1 TOPS mode data, covering the entire workflow from raw data selection, preprocessing, to final interferogram generation, phase unwrapping, and geocoding. The toolkit integrates the flexibility of Python scripts with the professional InSAR processing capabilities of GMTSAR, supporting parallel processing of multiple subswaths (IW1/IW2/IW3), batch interferogram pair generation, phase unwrapping, and result format conversion. It is suitable for both scientific research and engineering-level InSAR data processing.



# Installation
## Basic
- `Python` = 3.11.0
  -  $ wget https://mirrors.tuna.tsinghua.edu.cn/anaconda/miniconda/Miniconda3-py311_23.11.0-1-Linux-x86_64.sh
  -  $ bash Miniconda3-py311_23.11.0-1-Linux-x86_64.sh
- `gmt` = 6.4.0
- `gmtsar` = 6.6 (snaphu uses version 2.0.7)

## External Modules
- `numpy` = 1.26.4
- `matplotlib` = 3.10.8
- `scipy` = 1.16.3
- `beautifulsoup4` = 4.14.3
- `lxml` = 6.0.2
- `shapely` = 2.1.2
- `pandas` = 2.3.3
- `pyproj` = 3.7.2
- `geojson` = 3.2.0
- `gdal` = 3.12.1
- `libmambapy` = 2.5.0
- `netcdf4` = 1.7.3
- `xarray` = 2025.12.0

## Personal modules
- `pSAR`
- `pGMT`
- `pGMT5SAR`
- `pDATA`
- `pS1`
- `quadtree`
- `utm`

## python_script
pSAR_s1sorting.py\n
pSAR_S1select.py\n
pSAR_copdemdownload.py
pSAR_COP30M_downloader.py
pSAR_COP30M_tar2grd.py
pSAR_COP30M_untar.py
pSAR_imgformat.py
pSAR_unzip_parallel.py
pSAR_gmtsar_rawdir2prms.py
pSAR_gmtsar_raw2baseline.py
pSAR_gmtsar_dir2datalist.py
pSAR_gmtsar_tiff2slcs_paral.py
pSAR_gmtsar_dir2baseline.py
pSAR_baseline.py
pSAR_gmtsar_baseline2intfin.py
pSAR_gmtsar_refineSLC.py
gmtsar_unwrap.py
pSAR_gmtsar_dir2losvecs.py
pSAR_gmtsar_s1insar2roi.py
pSAR_gmtsar_dir2roi.py
pSAR_rasterio_fillnodata.py
pSAR_gmtsar_los2projvec.py
pSAR_orbitcor.py ##Orbit and terrain correction
pSAR_ui4poly.py ##Get coordinates
pSAR_netcdf4rsc_updateing.py ##Generate .rsc file
insar_preprocessing.py ##downsample

## cshell scripts required for interferometry
pop_config.csh
ginsar_unzip.sh
GMTSAR_s1_createTOPSframes.csh
gmtsar_preproc_batch_tops_esd.csh
gmtsar_intf_tops.csh
slc2amp_MUL.csh
gmtsar_filter.csh
gmtsar_geocode.csh
gmtsar_proj_ra2ll.csh

## scripts plot results
gmt_intf_plot.sh
fig_insat_gmt.sh ##Plot the LOS (GRD)file
gmt_plot_downsampled_polys.sh ##Plot the downsampled dataset
