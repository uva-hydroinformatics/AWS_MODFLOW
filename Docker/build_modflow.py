# -*- coding: utf-8 -*-
"""
Created on Fri Mar 31 17:19:25 2017

@author: wzell
"""

import os,subprocess,fiona,shutil,rasterio,itertools
import numpy as np
import flopy.modflow as mf

# --- SCRIPT PARAMETER SET START ---

# Tell the model where the MODFLOW binary/engine is. The current configuration assumes
# that the binary is in the MODFLOW directory.
binary_path = 'MODFLOW-NWT_64'

# Assign the model grid cell resolution (m). This model workflow assumes
# square model cells
dx_dy = 300

# Configure the cache for gdal subprocess calls. If you have GDAL installed
# (and in your path) you can build the MODFLOW framework arrays directly from
# the GIS layers.  If you do not have GDAL installed, set build_from_gis = False,
# in which case the model will be build from the ASCII files that I have already
# built in the Framework directory.
build_from_gis = True
default_cachemax = 4000
default_cache_config = ['--config','GDAL_CACHEMAX',str(default_cachemax)]

# --- SCRIPT PARAMETER SET STOP ----

# === HELPER FUNCTIONS AND CLASSES START =======

def write_raster_array(raster_in,file_out,fmt=None,multiplier=1,force_dims=None):
    '''Writes an individual raster to file.
    force_dims: (nrow,ncol) tuple to which ASCII array is reduced.'''
    
    with rasterio.open(raster_in,'r') as r:        
        iarray = r.read()[0,:,:] * multiplier
    
    if (force_dims is not None):
        force_rows,force_cols = force_dims
        irows,icols = np.shape(iarray)
        
        if (irows != force_rows) or (icols != force_cols):
            print '\tModel  dimensions = %i rows %i cols.' %(force_rows,force_cols)
            print '\tRaster dimensions = %i rows %i cols.' %(irows,icols)
            print 'Forcing raster to model dimensions.'
        
        if (irows > force_rows):
            iarray = iarray[0:force_dims[0],:]
            print 'Reduced array rows from %i to %i.' %(irows,np.shape(iarray)[0])
            
        if (icols > force_cols):
            iarray = iarray[:,0:force_dims[1]]
            print 'Reduced array rows from %i to %i.' %(icols,np.shape(iarray)[1])
        
        if (irows < force_rows):
            igap = int(force_rows) - int(irows)
            ifill = iarray[-1,:]
            for iadd in range(igap):
                iarray = np.vstack([iarray,ifill])
                
            print 'Added array rows; new nrow = %i.' %(np.shape(iarray)[0])
                
        if (icols < force_cols):
            jgap = int(force_cols) - int(icols)
            jfill = iarray[:,-1]
            for jadd in range(jgap):                
                iarray = np.column_stack((iarray,jfill))
                
            print 'Added array columns; new ncol = %i.' %(np.shape(iarray)[1])
                    
    np.savetxt(file_out,iarray,fmt=fmt)
        
    return 

def raster_to_model(data_fin=None,clipped_temp=None,raster_fout=None,\
                    bounds=None,delr=None,delc=None,\
                    model_epsg=None,resample='average',gdalwarp_path=None,\
                    cache_config=default_cache_config,cachemax=default_cachemax,\
                    save_temp=None,dstnodata=None,overwrite=True):
    '''Clips and resamples a raster to the model grid.'''

    if (gdalwarp_path is None):
        gdalwarp_path = 'gdalwarp'
        
    print 'Executing gdalwarp path: %s' %(gdalwarp_path)

    if (bounds is not None):
        
        # Clip the raster to the model domain bound
        print '\nClipping the raster to the model domain.\n'
        
        x_min,y_min,x_max,y_max = bounds
        
        clip_dem_cmd = [gdalwarp_path] + cache_config + ['-wm',str(cachemax/2)] + \
                       ['-te',str(x_min),str(y_min),str(x_max),str(y_max),data_fin,clipped_temp]
                       
        if overwrite:
            try:
                print '\tDeleting . . .'
                subprocess.call(['gdalmanage','delete',clipped_temp])
            except:
                print '\tDelete for overwrite failed.'
            
        subprocess.call(clip_dem_cmd)
    
        if (save_temp is not None):
            shutil.copyfile(clipped_temp,save_temp)
        
        resample_in = clipped_temp

    else:
        resample_in = data_fin        
        
    # Resample the raster to the model resolution
    print '\nResampling the raster to model grid resolution.\n'
    resample_cmd = [gdalwarp_path] + cache_config + ['-wm',str(cachemax/2)] + \
                   ['-r',resample,'-tr',str(delc),str(delr),\
                    resample_in,raster_fout]
        
    if dstnodata is not None:
        resample_cmd = resample_cmd + ['-dstnodata',str(dstnodata)]

    if overwrite:
        try:
            print '\tDeleting . . .'
            subprocess.call(['gdalmanage','delete',raster_fout])
        except:
            print '\tDelete for overwrite failed.'      
        
    subprocess.call(resample_cmd)

    return

def get_shp_extents(shp_fin=None):
    '''
    Returns the bounding box extents for a shapefile.
    '''

    with fiona.open(shp_fin,'r') as src:
        
        x_min,y_min,x_max,y_max = src.bounds
        
    return x_min,y_min,x_max,y_max

def shp_to_model(shp_fin=None,raster_fout=None,grid_dx=None,grid_dy=None,rasterize_field=None,\
                 no_data=-99999,cache_config=default_cache_config,burn_constant=None,dtype='Float64',\
                 bounds=None):
    '''Rasterizes a PROJECTED shapefile. Default: rasterize based upon user-specified shapefile
    field. If burn_constant is provided, burn that constant. If bounds are provided,
    clip the raster to those bounds.'''

    print '\nRasterizing shapefile: %s' %(shp_fin)
    print 'Writing output raster to: %s\n' %(raster_fout)

    layer_name = os.path.basename(shp_fin).replace('.shp','')
    rasterize_cmd = ['gdal_rasterize'] + cache_config

    if (burn_constant is not None):
        rasterize_cmd = rasterize_cmd + ['-a_nodata',str(no_data),'-burn',str(burn_constant),'-of','GTiff','-ot',dtype,'-l',layer_name,'-tr',str(grid_dx),str(grid_dy),shp_fin,raster_fout]
    else:            
        rasterize_cmd = rasterize_cmd + ['-a_nodata',str(no_data),'-a',rasterize_field,'-of','GTiff','-ot',dtype,'-l',layer_name,'-tr',str(grid_dx),str(grid_dy),shp_fin,raster_fout]

    if bounds is not None:
        rasterize_cmd = rasterize_cmd + ['-te'] + [str(x) for x in bounds]
    
    subprocess.call(rasterize_cmd)

    return

class Paths(object):
    '''
    This class aggregates and generates all of the file paths.
    
    wbd = Active model domain from Watershed Boundary Dataset
    dem = Elevation (cm) from National Elevation Dataset
    rch = Recharge (cm/year) from National Recharge Dataset (Reitz and Sanford)
    '''
    
    def __init__(self,model_name='James_Rivanna'):
        
        self.model_name = model_name        
        
        # Input data
        # ----------
        self.data_dir = 'Data'
        
        self.wbd_data     = os.path.join(self.data_dir,model_name + '_5070.shp')
        self.surfgeo_data = os.path.join(self.data_dir,'VA_SurficialGeology.shp')
        self.dem_data     = os.path.join(self.data_dir,'VA_NED.tif')
        self.rch_data     = os.path.join(self.data_dir,'VA_Recharge.tif')
        
        # Framework files. These are ASCII arrays derived from the GIS data.
        # MODFLOW reads these files in order to define the model grid,
        # system states (e.g., starting heads), and system properties
        # (e.g., hydraulic conductivity)
        # ---------------------------------------------------------------------
        self.model_frame_dir = 'Framework'
        
        # Intermediate rasters
        self.ibound_tif              = os.path.join(self.model_frame_dir,model_name + '_IBOUND.tif')
        self.model_dem_tif           = os.path.join(self.model_frame_dir,model_name + '_DEM.tif')
        self.model_rch_tif           = os.path.join(self.model_frame_dir,model_name + '_RCH.tif')
        
        self.landsurface_elev_file   = os.path.join(self.model_frame_dir,model_name + '.landsurface_elev')
        self.rch_file                = os.path.join(self.model_frame_dir,model_name + '.rch')
        self.surfgeo_file            = os.path.join(self.model_frame_dir,model_name + '.surfgeo')
        
        self.modeltop_file           = os.path.join(self.model_frame_dir,model_name + '.top')
        self.bottoms_file            = os.path.join(self.model_frame_dir,model_name + '.bottoms')
        self.ibound_file             = os.path.join(self.model_frame_dir,model_name + '.ibound')        
        self.starting_heads_file     = os.path.join(self.model_frame_dir,model_name + '.startheads')
                
        # MODFLOW input files (MODFLOW refers to these as 'packages')
        # -----------------------------------------------------------
        self.modflow_dir = 'MODFLOW'
        self.mf_version  = 'mfnwt'
        self.mf_bat_file = os.path.join(self.modflow_dir,model_name + '.' + self.mf_version + '.bat')
        
        self.nam_file   = os.path.join(self.modflow_dir,model_name + '.nam')
        self.upw_pkg    = os.path.join(self.modflow_dir,model_name + '.upw')
        self.rch_pkg    = os.path.join(self.modflow_dir,model_name + '.rch')
        self.drn_pkg    = os.path.join(self.modflow_dir,model_name + '.drn')
        self.oc_pkg     = os.path.join(self.modflow_dir,model_name + '.oc')
        
        # Scratch workspace
        # -----------------
        self.scratch_dir = 'Scratch'        
        
        self.dem_clipped = os.path.join(self.scratch_dir,model_name + '_temp_dem.tif')
        self.rch_clipped = os.path.join(self.scratch_dir,model_name + '_temp_rch.tif')
        
        # Make the MODFLOW and Framework directories if they don't exist.
        # Note that the data directory MUST already exist.
        for idir in [self.modflow_dir,self.model_frame_dir,self.scratch_dir]:
            if not os.path.exists(idir):
                os.makedirs(idir)
            
        return
        
class Frame(object):
    '''
    This class defines the model framework and includes the methods
    that generate the model input arrays from the GIS data.
    '''
    
    def __init__(self,Paths=None,dx_dy=None,model_epsg=5070,\
                 lay_thick=[100.],laytyp=[0],\
                 hdry=-888.,hnoflo=-999.):
        
        self.model_name = Paths.model_name
        self.model_epsg = model_epsg
        
        # Attributes that will be used by the MODFLOW .dis package
        # (as well as a few convenient derivatives)
        self.delr,self.delc = dx_dy,dx_dy
        self.cell_area = self.delr * self.delc        
        self.nlay = len(lay_thick)        
        self.lay_thick = lay_thick
        self.laytyp = laytyp                
        
        # Attributes that will be used by the MODFLOW .upw package
        self.hnoflo = hnoflo
        self.hdry = hdry
        
        # If the geospatial data has already been processed, grab
        # the arrays
        try:
            self.ibound = np.genfromtxt(Paths.ibound3D_file)
            self.landsurface = np.genfromtxt(Paths.landsurface_elev_file)
            self.nrow,self.ncol = np.shape(self.top)
            self.shape = (self.nlay,self.nrow,self.ncol)
        except:
            pass
        
        return
        
    def build_frame(self,Paths=None,\
                    rch_unit_multiplier=1./365.,dem_unit_multiplier=0.01):
        '''
        Writes model domain, IBOUND zones, and framework arrays for a regional
        model from national GIS data.
        
        Default conversion factors:
            default_rch_unit_multiplier = 1./365. (Convert m/year to m/day)
            default_dem_unit_multiplier = 0.01    (Convert cm to m) 
        '''

        self.bbox = get_shp_extents(shp_fin=Paths.wbd_data)
        
        # GENERATE MODEL RESOLUTION RASTERS
        # ---------------------------------------
        
        # Generate the IBOUND array
        shp_to_model(shp_fin=Paths.wbd_data,raster_fout=Paths.ibound_tif,\
                     grid_dx=self.delc,grid_dy=self.delr,bounds=self.bbox,\
                     no_data=0,burn_constant=1,dtype='Int32')
                     
        # Clip the DEM and recharge information to the model grid
        raster_to_model(data_fin=Paths.dem_data,clipped_temp=Paths.dem_clipped,\
                        raster_fout=Paths.model_dem_tif,bounds=self.bbox,resample='average',\
                        delr=self.delr,delc=self.delc,model_epsg=self.model_epsg,\
                        dstnodata=-9999.)
                        
        raster_to_model(data_fin=Paths.rch_data,clipped_temp=Paths.rch_clipped,\
                        raster_fout=Paths.model_rch_tif,bounds=self.bbox,resample='average',\
                        delr=self.delr,delc=self.delc,model_epsg=self.model_epsg,\
                        dstnodata=-9999.)

        # GENERATE THE ASCII ARRAYS. Note that the dimensions of the IBOUND
        # array are used as the template for all other arrays and used to
        # define (nrow,ncol) for MODFLOW
        # -------------------------------
        write_raster_array(Paths.ibound_tif,Paths.ibound_file,fmt='%12i')
        self.nrow,self.ncol = np.shape(np.genfromtxt(Paths.ibound_file))

        write_raster_array(Paths.model_dem_tif,Paths.landsurface_elev_file,fmt='%15.8e',\
                           multiplier=dem_unit_multiplier,force_dims=(self.nrow,self.ncol))                                      
                        
        write_raster_array(Paths.model_rch_tif,Paths.rch_file,fmt='%15.8e',\
                           multiplier=rch_unit_multiplier,force_dims=(self.nrow,self.ncol))
        
        # Now capture the ASCII files as numpy arrays that can be passed to the 
        # Flopy model constructors
        self.ibound = np.genfromtxt(Paths.ibound_file)
        self.top    = np.genfromtxt(Paths.landsurface_elev_file)
        self.bottom = self.top - self.lay_thick
        self.rch    = np.genfromtxt(Paths.rch_file)
                        
        return

def build_drain_input(mfFrame=None,stages=None,condmult=1):
    '''
    Generates a dictionary of drain specifiers for MODFLOW input.
    Each drain in each stress period must be defined with 
    [lay, row, col, stage, cond]
    '''

    conductance = mfFrame.cell_area * condmult

    drns = []    
    for irow,icol in itertools.product(range(mfFrame.nrow),range(mfFrame.ncol)):
        
        if np.isfinite(stages[irow,icol]):        
            drns.append([0,irow,icol,stages[irow,icol],conductance])
    
    return {0:drns}

# === HELPER FUNCTIONS AND CLASSES STOP ========

def main():
    '''
    This is the main function.
    '''

    # Package all the required file paths into the Paths object        
    mfPaths = Paths()
    
    # Package all the required framework specifications into the mfFrame object
    mfFrame = Frame(Paths=mfPaths,dx_dy=dx_dy)
    
    if build_from_gis:
        
        # Build the model framework ASCII files from the GIS layers. Note that this
        # requires a GDAL installation.  If you don't want to get into that you
        # can skip this step and simply build the model from the ASCII files that I've
        # already created.
        mfFrame.build_frame(Paths=mfPaths)
    # ---------------------------------------------
    # ---------------------------------------------
    # Now use Flopy to build the MODFLOW model packages
    # ---------------------------------------------
    # ---------------------------------------------
    
    start_dir = os.getcwd()
    os.chdir(mfPaths.modflow_dir) # This is simplest if done inside the MODFLOW directory
    
    # Initialize a Flopy model object. This is the base class around which the model
    # packages are built.
    Modflow = mf.Modflow(mfFrame.model_name,external_path='./',version=mfPaths.mf_version)
    
    # The .oc ('output control') package specifies how the model output is written.
    # This model includes a single steady state stress period. Save the
    # distribution of heads as well as the flow budget/mass balance to binaries.
    # These can be plotted or converted to rasters (the current version of the script
    # doesn't do any post-processing.)
    oc = mf.ModflowOc(Modflow,stress_period_data={(0,0):['SAVE HEAD','SAVE BUDGET']})
    
    # The .dis and .bas packages define the model framework. I've already defined
    # the framework attributes using the mfFrame object and simply pass those
    # attributes to the constructor.
    dis = mf.ModflowDis(Modflow,mfFrame.nlay,mfFrame.nrow,mfFrame.ncol,\
                         delr=mfFrame.delr,delc=mfFrame.delc,\
                         top=mfFrame.top,botm=mfFrame.bottom)
                         
    bas = mf.ModflowBas(Modflow,ibound=mfFrame.ibound,strt=mfFrame.top,hnoflo=mfFrame.hnoflo)
    
    # The .upw package describes the system properties (e.g., transmissivity/conductivity).
    # For this model I simply give it a constant hydraulic conductivity field. This model
    # converges but I have no idea how physically realistic it is. If you would
    # like to make it more physically realistic (e.g., try to fit head or discharge
    # data) you would need to estimate the hydraulic conductivity field via
    # calibration/inverse modeling
    hk = np.ones(np.shape(mfFrame.ibound))
    upw = mf.ModflowUpw(Modflow,laytyp=mfFrame.laytyp,hk=hk)
    
    # The .nwt package defines the solver specs. Just use the defaults.
    nwt = mf.ModflowNwt(Modflow)
    
    # RECHARGE INPUTS TO THE SYSTEM
    # -----------------------------
    # The .rch packages specifies recharge/precipitation inputs to the water table.
    # Remember that I have already generated an array from the GIS layer and attached
    # it to the mfFrame object.
    rch = mf.ModflowRch(Modflow,nrchop=3,rech={0:mfFrame.rch})
    
    # BASEFLOW DISCHARGE FROM THE SYSTEM
    # ----------------------------------
    # The .drn package is one method of simulating the discharge of groundwater as
    # base-flow in streams in rivers.  Define every landsurface cell as a drain
    # in order to allow the discharge network to emerge from topography.
    drn_stages = mfFrame.top
    drn_stages[mfFrame.ibound.squeeze() <= 0] = np.nan
    drn_input = build_drain_input(mfFrame=mfFrame,stages=drn_stages)
    
    drn = mf.ModflowDrn(Modflow,stress_period_data=drn_input)
    
    # Now write the files.  Flopy can also run the model if you tell it where the
    # binary is, but if I understood your method correctly you will be invoking something
    # from hydroshare. For convenience I am writing a windows .bat file that
    # can be used to run the model.
    Modflow.write_input()    
    os.chdir(start_dir)
    
    with open(mfPaths.mf_bat_file,'w') as fout:
        fout.write('%s %s' %(binary_path,os.path.basename(mfPaths.nam_file)))

    return       
main()
