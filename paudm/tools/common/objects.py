
"""
PAUdm objects definition

@author: s.serrano_

To constrain the use of memory to one CCD (Image), process mosaics using the following structure:


in_mosaic = objects.Mosaic(path='/Users/santi/', name='in.fits', mode='open')
out_mosaic = objects.Mosaic(path='/Users/santi/', name='out.fits', mode='new', global_header=in_mosaic.header)

for ccd_id in range(1, in_mosaic.get_keyword('NEXTEND')+1):
    current_ccd = in_mosaic.load_ccd(ccd_id)
    out_mosaic.append_ccd(current_ccd)


"""

import logging
log = logging.getLogger('paudm.tools.common.objects')
import os
import datetime
import time
import re
import math
import yaml
import numpy as np
import pyfits
from scipy.interpolate import interp1d
from astropysics.coords.coordsys import AngularCoordinate
from paudm.tools.common import constants
from pkg_resources import resource_listdir, resource_string, resource_filename
config = {}
from paudm.tools.db import model

#instrument = yaml.safe_load(resource_string('paudm.resources.instrument.pau','general.yaml'))




class Mosaic(object):
    '''
    Defines a standard Mosaic object.
    
    @param: path, name, mode (open/new) [global_header]
    
    Attributes:
      - name
      - path
      - header
      - assoc_mosaic (raw, mbias, mflat, mask, weight, cosm,..)
      - assoc_db
      
    Built-in Methods:
      - value = get_keyword(key)
      - status = set_keyword(key, value)
      - image = load_image(image_id)
      - status = store_image(image)
      - status = append_image(image_id)
      - status = store_global_header()
      - status = db_insert( )
    '''
    
    
    def __init__(self, path, filename, mode, global_header = None):
        # Mosaic Constructor
        self.filename = filename
        self.path = path
        self.header = global_header # TODO: Turn into a property
        self.assoc_mosaic = {}
        self.assoc_db = None
        self.badcols = {}
        self.pyfits_mosaic = None
        
        if mode.lower() == 'open':
            # Load global header
            self.header = pyfits.getheader(os.path.join(self.path, self.filename))
            self.pyfits_mosaic = pyfits.open(os.path.join(self.path, self.filename))
        elif mode.lower() == 'new':
            # Remove if exists
            if os.path.exists(os.path.join(self.path, self.filename)):
                os.remove(os.path.join(self.path, self.filename))
            # Create new Primary fits
            new_hdu = pyfits.PrimaryHDU()
            # Use blank header if not given
            if self.header == None:
                self.header = new_hdu.header
            pyfits.writeto(os.path.join(self.path, self.filename), data = np.array([], dtype = 'uint8'), header = self.header)
        
    
    # TODO: Turn into the getter of self.header
    def get_keyword(self, key):
        # Get keyword from global header
        if not self.header.has_key(key):
            raise KeyError('Keyword does not exist in header')
        else:
            return self.header.get(key)
    
    # TODO: Turn into the setter of self.header
    def set_keyword(self, key, value, comments = None):
        # Set keyword in global header
        if comments == None:
            self.header.update(key, value)
        else:
            self.header.update(key, value, comments)
    
    # TODO: Turn into the deletter of self.header
    def remove_keyword(self, key):
        # Removes keyword from header
        del self.header[key]
    
    def load_image(self, image_id):
        # Build and return image from Mosaic
        return Image(mode = 'open', path = self.path, filename = self.filename, extension = image_id, parent_mosaic = self)
    
    def store_image(self, image):
        # Overrite an image in Mosaic
        pyfits.update(str(os.path.join(self.path, self.filename)), data = image.data, header = image.header, ext = image.id)
    
    def append_image(self, image):
        # Append an image to Mosaic
        pyfits.append(str(os.path.join(self.path, self.filename)), data = image.data, header = image.header)
    
    def append_image_fast(self, image):
        # Manual appending with regular i/o files (faster than Pyfits)
        hdu_tmp = pyfits.ImageHDU(image.data)
        hdu_tmp.header = image.header
        hdu_tmp.writeto(os.path.join(self.path, self.filename) + '.tmp', clobber = True)
        hdu_fid = open(os.path.join(self.path, self.filename) + '.tmp','r')
        hdu_raw = hdu_fid.read()
        # Find beggining of XTENSION section (Image HDU)
        xtension_ipos = hdu_raw.find('XTENSION')
        hdu_fid.close()
        # Append to Mosaic file
        hdu_mosaic = open(os.path.join(self.path, self.filename),'a')
        hdu_mosaic.write(hdu_raw[xtension_ipos : ])
        hdu_mosaic.close()
        # Clean
        os.remove(os.path.join(self.path, self.filename) + '.tmp')
        del hdu_raw
        del hdu_tmp  
    
    def store_global_header(self):
        # Save changes in global header
        pyfits.update(os.path.join(self.path, self.filename), data = np.array([], dtype = 'uint8'), header = self.header)
    
    # DATA BASE methods
    def fields_from_header(self):
        # Get mosaic header fields needed in the data base
        header = self.header
        fields = {}
        if 'KIND' in header:
            fields['kind']       = header['KIND']
        else:
            fields['kind']       = header['OBSTYPE']
        fields['exp_num']    = header['EXPNUM']
        fields['obs_title']  = header['OBJECT']
        fields['ra']         = AngularCoordinate(header['RA'], sghms=True).degrees
        fields['dec']        = AngularCoordinate(header['DEC'], sghms=False).degrees
        fields['equinox']    = header['EQUINOX']
        fields['date_obs']   = datetime.date(*(time.strptime(header['DATE-OBS'],'%Y-%m-%d')[:3]))
        fields['time_obs']   = datetime.time(*(time.strptime(header['UTSTART'].split('.')[0],'%H:%M:%S')[3:6]))
        fields['date_creat'] = datetime.date(*(time.strptime(header['DATE-OBS'],'%Y-%m-%d')[:3]))
        fields['time_creat'] = datetime.time(*(time.strptime(header['UTSTART'].split('.')[0],'%H:%M:%S')[3:6]))
        fields['rjd_obs']    = header['RJD-OBS']
        fields['exp_time']   = header['EXPTIME']
        fields['airmass']    = header['AIRMASS']
        fields['telfocus']   = header['TELFOCUS']
        fields['instrument'] = header['INSTRUME']
        if header['FTRAY'] != 'None':
            fields['filtertray'] = header['FTRAY']
        else:
            fields['filtertray'] = None
        if header['FTRAYTMP'] != 'None':
            fields['filtertray_tmp'] = header['FTRAYTMP'] 
        else:
            fields['filtertray_tmp'] = None
        fields['guide_enabled']  = (header['AGUIDER'] == 'T')
        fields['guide_fwhm']     = header['GSFWHM']
        fields['guide_var']      = header['GSVAR']
        
        fields['nextend']  = header['NEXTEND']
        fields['wind_spd'] = header['WSPEED']
        fields['wind_dir'] = header['WDIRECT']
        fields['amb_temp'] = header['TEMPOUT']
        fields['humidity'] = header['HUMIDITY']
        fields['pressure'] = header['PRESSURE']
        
        fields['eqa_1'] = header['EQA01']
        fields['eqa_2'] = header['EQA02']
        fields['eqa_3'] = header['EQA03']
        fields['eqa_4'] = header['EQA04']
        fields['eqa_5'] = header['EQA05']
        
        return fields
    
    # NIGHTLY Delegates
    def create_overscanned_mosaic(self):
        import paudm.pipeline.nightly.delegates
        return paudm.pipeline.nightly.delegates.mosaic_create_overscanned_mosaic(self)
    
    def load_ccd_raw(self, instrument, ccd_num, correct_overscan = True, max_overscan_dev = None):
        # Load amps and build overscan-corrected CCD
        import paudm.pipeline.nightly.delegates
        return paudm.pipeline.nightly.delegates.load_ccd_raw(self, instrument, ccd_num, correct_overscan, max_overscan_dev)
    
    # PIXELSIM Delegates
    def initialize_sim_params(self, exposure, environment, config, instrument):
        import paudm.pipeline.pixelsim.delegates
        return paudm.pipeline.pixelsim.delegates.initialize_sim_params(self, exposure, environment, config, instrument)
    
    def simulate(self,config, instrument ):
        import paudm.pipeline.pixelsim.delegates
        return paudm.pipeline.pixelsim.delegates.simulate_mosaic(self, config, instrument)
    
    def initialize_global_header(self, config, AMPS_X_CCD):
        import paudm.pipeline.pixelsim.delegates
        return paudm.pipeline.pixelsim.delegates.initialize_global_header(self, config['common']['general']['project'].lower(), config['release'], config['ccd_limit'], AMPS_X_CCD, config['post_prod'])
    
    def infer_pixel_rectangle(self, config, instrument):
        import paudm.pipeline.pixelsim.delegates
        return paudm.pipeline.pixelsim.delegates.infer_pixel_rectangle(self, config, instrument)
    


class Image(object):
    '''
    Defines a standard Image object.
    
    @param: path, name, extension, mode (open/new) [image_header]
    
    Attributes:
      - id
      - header
      - data[]
      
    Built-in Methods:
      - value = get_keyword(key)
      - status = set_keyword(key, value)
      - status = db_insert( )
    '''
    
    def __init__(self, mode, path = None, filename = None, extension = None, image_header = None, parent_mosaic = None):
        # Image Constructor
        self.id = extension
        self.header = image_header # TODO: Turn into a property
        self.data = None
        self.amp_array = None
        self.parent_mosaic = parent_mosaic
        
        # ORM object stored in the database
        self._db_object = None
        
        if mode == 'open':
            # Load image
            self.py_image = self.parent_mosaic.pyfits_mosaic[self.id]
            self.header = self.py_image.header
            self.data = self.py_image.data
        elif mode == 'new':
            # New image
            # Create new Primary fits
            new_hdu = pyfits.ImageHDU()
            self.data = new_hdu.data
            # Use blank header if not given
            if self.header == None:
                self.header = new_hdu.header
    
    # TODO: Turn into the getter of self.header
    def get_keyword(self, key):
        # Get keyword from image header
        if not key in self.header.keys():
            raise KeyError('Keyword does not exist in header')
        else:
            return self.header.get(key)
    
    # TODO: Turn into the setter of self.header
    def set_keyword(self, key, value, comments = None):
        # Set keyword in global header
        if comments == None:
            self.header.update(key, value)
        else:
            self.header.update(key, value, comments)
    
    # TODO: Turn into the deletter of self.header
    def remove_keyword(self, key):
        # Removes keyword from header
        del self.header[key]
    
    def save_check(self, filename):
        # Save a single extension file with image data for visual check
        if os.path.exists(filename):
            os.remove(filename)
        pyfits.writeto(filename, data = self.data, header = self.header)
    
    
    # DATA BASE Methods
    def fields_from_header(self, pixel_scale):
        fields = {}
        #mapping
        fields['image_num']  = self.header['IMAGEID']
        fields['ccd_num']    = self.header['CCDNUM']
        fields['amp_num']    = self.header['AMPNUM']
        if self.header['FILTER'] != 'None':
            fields['filter'] = self.header['FILTER']
        else:
            fields['filter'] = None
        fields['wavelength'] = self.header['WAVELEN']
        fields['waveband']   = self.header['WAVEBAND']
        fields['gain']       = self.header['GAIN']
        fields['rdnoise']    = self.header['RDNOISE']
        fields['naxis1']     = self.header['NAXIS1']
        fields['naxis2']     = self.header['NAXIS2']
        fields['cqa_1']      = self.header['CQA01']
        fields['cqa_2']      = self.header['CQA02']
        fields['cqa_3']      = self.header['CQA03']
        fields['cqa_4']      = self.header['CQA04']
        fields['cqa_5']      = self.header['CQA05']
        
        # Calculate Sky Corners of the image
        if self.parent_mosaic.header['OBSTYPE'] in ['TARGET', 'RED_SCI', 'RED_MASK', 'RED_WEIGHT']:
            from paudm.pipeline.pixelsim import wcsUtils
            from paudm.pipeline.pixelsim import simUtils
            
            # Load CRVAL
            image_wcs = wcsUtils.wcsType(sky_ref=simUtils.skyPoint(ra=self.header['CRVAL1'], dec=self.header['CRVAL2']),
                                         pix_ref=simUtils.pixelPoint(0, 0), # Dummy Value
                                         rot_degrees=0,
                                         instrument_PIXEL_SCALE = pixel_scale) # Dummy Value
            # Override WCS
            image_wcs.CRPIX1 = self.header['CRPIX1']
            image_wcs.CRPIX2 = self.header['CRPIX2']
            image_wcs.CD1_1 = self.header['CD1_1']
            image_wcs.CD1_2 = self.header['CD1_2']
            image_wcs.CD2_1 = self.header['CD2_1']
            image_wcs.CD2_2 = self.header['CD2_2']
            
            ccd_corners = [image_wcs.pix2sky(simUtils.pixelPoint(0, 0)),
                           image_wcs.pix2sky(simUtils.pixelPoint(self.header['NAXIS1'], 0)),
                           image_wcs.pix2sky(simUtils.pixelPoint(0, self.header['NAXIS2'])),
                           image_wcs.pix2sky(simUtils.pixelPoint(self.header['NAXIS1'], self.header['NAXIS2']))]
            
            fields['ra_min']    = min([element.ra for element in ccd_corners])
            fields['ra_max']    = max([element.ra for element in ccd_corners])
            fields['dec_min']   = min([element.dec for element in ccd_corners])
            fields['dec_max']   = max([element.dec for element in ccd_corners])
            
        return fields
    
    
    # NIGHTLY Methods
    def bias_subtract(self):
        import paudm.pipeline.nightly.delegates
        return paudm.pipeline.nightly.delegates.bias_subtract(self)
    
    def flatfield_correct(self):
        import paudm.pipeline.nightly.delegates
        return paudm.pipeline.nightly.delegates.flatfield_correct(self)
    
    def badcols_interpolate(self):
        import paudm.pipeline.nightly.delegates
        return paudm.pipeline.nightly.delegates.badcols_interpolate(self)
    
    def badpix_interpolate(self):
        import paudm.pipeline.nightly.delegates
        return paudm.pipeline.nightly.delegates.badpix_interpolate(self)
    
    # PIXELSIM Methods
    def initialize_extension_header(self, config, instrument):
        import paudm.pipeline.pixelsim.delegates
        return paudm.pipeline.pixelsim.delegates.initialize_extension_header(self, config, instrument)
    
    def update_simulation_keywords(self):
        import paudm.pipeline.pixelsim.delegates
        return paudm.pipeline.pixelsim.delegates.update_simulation_keywords(self)
    
    def set_sim_WCS(self, telescope_pointing, config,instrument):
        import paudm.pipeline.pixelsim.delegates
        return paudm.pipeline.pixelsim.delegates.set_sim_WCS(self, telescope_pointing, config, instrument)
    
    def contains_pixels( self, pixelPos, buffer_width=0 ):
        import paudm.pipeline.pixelsim.delegates
        return paudm.pipeline.pixelsim.delegates.contains_pixels( self, pixelPos, buffer_width=0 )
    
    def prepare_psf(self):
        import paudm.pipeline.pixelsim.delegates
        return paudm.pipeline.pixelsim.delegates.prepare_psf(self)
    
    def simulate_image(self, config, instrument):
        import paudm.pipeline.pixelsim.delegates
        return paudm.pipeline.pixelsim.delegates.simulate_image(self, config, instrument)
    
    def post_production(self, config, instrument):
        import paudm.pipeline.pixelsim.delegates
        return paudm.pipeline.pixelsim.delegates.post_production(self, config, instrument)
    
    def get_sim_amp(self, amp_num, config, instrument):
        import paudm.pipeline.pixelsim.delegates
        return paudm.pipeline.pixelsim.delegates.get_sim_amp(self, amp_num, config, instrument)
    


class Catalogue(object):
    '''
    Defines a standard Catalogue object.
    
    @param: TBD
    
    Attributes:
      - path
      - filename
      - objects[]
      
      usage example: 
        objects[obj_num] -> list of all object properties for obj_num
        objects[obj_num][property] -> value of property for obj_num
        objects[property] -> list property value for all objects
        objects._names -> list all properties 
    
    Built-in Methods:
      - status = db_ingest( )
    '''
    
    def __init__(self, objs = None, ext = None, p = None, f = None):
        if objs is None:
            self.objects = []
        else:
            self.objects = objs
        self.query = None
        self.path = p
        self.filename = f
        self.extension = ext
        
    
    def load_all_query(self):
        # Load all objects of the query to the objects list
        self.objects = self.query.all()
    
    def sync_with_db(self, session):
        # Synchronize the objects list with the data base
        if len(self.objects) == 0:
            log.warning("No objects in the catalogue list. Try load_all_query?")
        for db_element in self.objects:
            if db_element.sed_type == None or db_element < 0:
                log.warning("Why this is here???")
                log.warning("%s"%vars(db_element))
                print "-------!!!!!!!!!!!!!!!!------------"

            session.add(db_element)
            session.flush()
    
    def merge_with_db(self, session):
        # Synchronize the objects list with the data base
        if len(self.objects) == 0:
            log.warning("No objects in the catalogue list. Try load_all_query?")
        for db_element in self.objects:
            session.merge(db_element)
    
    # DETECTION TABLE methods
    @classmethod
    def create_db_detections_catalogue(cls, path, filename, extension = None, total_extensions = 1, input_type = 'ldac', mosaic = None):
        # Create a list of detections with the entries of the catalogue
        catalog = cls()
        catalog.path = path
        catalog.filename = filename
        total_extensions = total_extensions
        catalog.objects = list()
        catalog.mosaic = mosaic
        
        if input_type == 'ldac':
            # LDAC is regular SExtractor fits file for catalogues
            log.debug("Open_LDAC file and create a catalog Object of detections")
            n_extensions = 0
            while n_extensions < total_extensions:
                catalog.extension = extension #First CCD to compute
                detections = pyfits.getdata(os.path.join(catalog.path, catalog.filename), ext = 2 * catalog.extension)
                catalog._names = detections.dtype.names
                extension +=1
                n_extensions +=1
                
                for detection in detections:
                    dict_tmp = {}
                    for property_id in range(0, len(catalog._names)):
                        property = catalog._names[property_id]
                        dict_tmp[property] = detection[property_id]
                    fields = catalog.detections_mapping(dict = dict_tmp)
                    db_image = model.session.query(model.Image).filter_by(mosaic = mosaic.assoc_db, image_num = catalog.extension).one()
                    det_tmp = model.Detection(band = db_image.filter,
                                              production = mosaic.assoc_db.production,
                                              image = db_image,
                                              ccd = catalog.extension,
                                              zp_offset = 0.0,
                                              **fields)
                    catalog.objects.append(det_tmp)
        else:
            error_msg = "Invalid input_type: %s" %input_type
            log.error(error_msg)
            raise Exception, error_msg
        
        return catalog
    
    
    def detections_mapping(self, dict):
        # Map data base detections elements from SExtractor output
        
        fields = {}
        #mapping
        fields['background'] = float(dict['BACKGROUND'])
        fields['class_star'] = float(dict['CLASS_STAR'])
        fields['flux_auto'] = float(dict['FLUX_AUTO'])
        fields['flux_err_auto'] = float(dict['FLUXERR_AUTO'])
        fields['flux_psf'] = float(dict['FLUX_PSF'])
        fields['flux_err_psf'] = float(dict['FLUXERR_PSF'])
        fields['flux_model'] = float(dict['FLUX_MODEL'])
        fields['flux_err_model'] = float(dict['FLUXERR_MODEL'])
        fields['flags'] = int(dict['FLAGS']) + int(dict['IMAFLAGS_ISO'])
        fields['elongation'] = float(dict['ELONGATION'])
        fields['dec'] = float(dict['DELTAWIN_J2000'])
        fields['ra'] = float(dict['ALPHAWIN_J2000'])
        fields['x'] = float(dict['XWIN_IMAGE'])
        fields['y'] = float(dict['YWIN_IMAGE'])
        
        return fields
    
    
    # TRUTH_OBJECT TABLE methods
    @classmethod
    def create_db_truth_object_catalogue(cls, input_type = 'mice_galaxy_cat', truth_objects_list = [], production = None):
        # Create a list of truth objects
        """
        Available input types:
         - mice_galaxy_cat
         - star_cat        
        """
        catalog = cls()
        
        if input_type == 'mice_galaxy_cat':
            
            for gal in truth_objects_list:
                # Mapping between galaxyType and truth_object (galaxy)
                fields = catalog.sim_gal_mapping(gal)
                # Create truth_object db element
                gal_tmp = model.Truth_Object( stargalaxy = False,
                                              production = production,
                                              **fields)
                # Append to catalog.objects list
                catalog.objects.append(gal_tmp)
            
        elif input_type == 'star_cat':
            
            for star in truth_objects_list:
                
                # Mapping between starType and truth_object (star)
                fields = catalog.sim_star_mapping(star)
                # Create truth_object db element
                star_tmp = model.Truth_Object( stargalaxy = True,
                                               production = production,
                                               **fields)
                # Append to catalog.objects list
                catalog.objects.append(star_tmp)
                
        else:
            error_msg = "Invalid input_type: %s" %input_type
            log.error(error_msg)
            raise Exception, error_msg
        
        return catalog
        
    
    
    def sim_star_mapping(self, star):
        # Map data base truth_object elements from pixelsim starType
        
        fields = {}
        #mapping
        fields['ra'] = star.ra
        fields['dec'] = star.dec
        fields['mag_u'] = star.mags['u']
        fields['mag_g'] = star.mags['g']
        fields['mag_r'] = star.mags['r']
        fields['mag_i'] = star.mags['i']
        fields['mag_z'] = star.mags['z']
        fields['mag_y'] = star.mags['y']
        fields['mag_n01'] = star.mags['n01']
        fields['mag_n02'] = star.mags['n02']
        fields['mag_n03'] = star.mags['n03']
        fields['mag_n04'] = star.mags['n04']
        fields['mag_n05'] = star.mags['n05']
        fields['mag_n06'] = star.mags['n06']
        fields['mag_n07'] = star.mags['n07']
        fields['mag_n08'] = star.mags['n08']
        fields['mag_n09'] = star.mags['n09']
        fields['mag_n10'] = star.mags['n10']
        fields['mag_n11'] = star.mags['n11']
        fields['mag_n12'] = star.mags['n12']
        fields['mag_n13'] = star.mags['n13']
        fields['mag_n14'] = star.mags['n14']
        fields['mag_n15'] = star.mags['n15']
        fields['mag_n16'] = star.mags['n16']
        fields['mag_n17'] = star.mags['n17']
        fields['mag_n18'] = star.mags['n18']
        fields['mag_n19'] = star.mags['n19']
        fields['mag_n20'] = star.mags['n20']
        fields['mag_n21'] = star.mags['n21']
        fields['mag_n22'] = star.mags['n22']
        fields['mag_n23'] = star.mags['n23']
        fields['mag_n24'] = star.mags['n24']
        fields['mag_n25'] = star.mags['n25']
        fields['mag_n26'] = star.mags['n26']
        fields['mag_n27'] = star.mags['n27']
        fields['mag_n28'] = star.mags['n28']
        fields['mag_n29'] = star.mags['n29']
        fields['mag_n30'] = star.mags['n30']
        fields['mag_n31'] = star.mags['n31']
        fields['mag_n32'] = star.mags['n32']
        fields['mag_n33'] = star.mags['n33']
        fields['mag_n34'] = star.mags['n34']
        fields['mag_n35'] = star.mags['n35']
        fields['mag_n36'] = star.mags['n36']
        fields['mag_n37'] = star.mags['n37']
        fields['mag_n38'] = star.mags['n38']
        fields['mag_n39'] = star.mags['n39']
        fields['mag_n40'] = star.mags['n40']
        fields['mag_n41'] = star.mags['n41']
        fields['mag_n42'] = star.mags['n42']
        
        #fields['sed_type'] = star.sed
        fields['sed_em_line'] = -1
        #fields['redshift'] = 0.0 #redshift is nullable, so this is optional
        
        return fields
        
    
    def sim_gal_mapping(self, gal):
        # Map data base truth_object elements from pixelsim galType
        
        # A lot of the work is already done in the starType mapping 
        # since the starType has a subset of the galType fields
        fields = self.sim_star_mapping( gal )
        
        # and now add the galaxy characteristics
        fields['bulge2flux_ratio'] = gal.bulge2flux_ratio
        fields['bulge_length'] = gal.bulge_length
        fields['bulge_axis_ratio'] = gal.bulge_axis_ratio
        fields['bulge_angle'] = gal.bulge_position_angle
        fields['disk_length'] = gal.disk_length
        fields['disk_axis_ratio'] = gal.disk_axis_ratio
        fields['disk_angle'] = gal.disk_position_angle
        fields['redshift'] = gal.redshift
        
        return fields
        
    
    
    # Sync local db with usno stars db at sky_region
    def usno_load_external(self, sky_region, config):
    
        log.debug("Query to USNO external data base...")
        import urllib
        #url='http://irsa.ipac.caltech.edu/cgi-bin/Gator/nph-query'
        url = config['urls']['usno'][0]
        fmt='csv'
        
        dec_mid = (sky_region['dec']['min'] + sky_region['dec']['max'])/2.0
        
        separator = "+"
        
        sky_regions = []
        if sky_region['ra']['max'] > sky_region['ra']['min']:
            sky_regions.append(sky_region)
        else:
            sky_region1 = dict(ra=dict(min=sky_region['ra']['min'],max=360.0),dec=dict(min=sky_region['dec']['min'],max=sky_region['dec']['max']))
            sky_region2 = dict(ra=dict(min=0.0,max=sky_region['ra']['max']),dec=dict(min=sky_region['dec']['min'],max=sky_region['dec']['max']))
            sky_regions.append(sky_region1)
            sky_regions.append(sky_region2)
            
        
        for sky_reg in sky_regions:
            ra_mid = (sky_reg['ra']['min'] + sky_reg['ra']['max'])/2.0
            
            cmd = ("outfmt=1&objstr=%f%s%f&spatial=Polygon&polygon=%f%s%f,%f%s%f,%f%s%f,%f%s%f&selcols=ra,dec,e_ra,e_dec,b1_mag,b2_mag,r1_mag,r2_mag,USNO_B1&catalog=usno_b1" %(ra_mid, separator, dec_mid, sky_reg['ra']['min'],separator,sky_reg['dec']['min'],sky_reg['ra']['max'],separator,sky_reg['dec']['min'],sky_reg['ra']['max'],separator,sky_reg['dec']['max'],sky_reg['ra']['min'],separator,sky_reg['dec']['max']))
            #log.debug( cmd );
            
            usno_out = (urllib.urlopen(url+'?%s' % cmd)).read()
            
            log.debug("...selection received.")
            
            #print "USNO OUTPUT:"
            #print usno_out
            usno_catalogue = usno_out.split('\n')
            re_comm =  re.compile("\\\\|\|") ##if( re.match("\\\\", usno_line[0]) is not None  or re.match("\|", usno_line[0]) is not None )
            re_html = re.compile("html")
            for usno_line in usno_catalogue:
                if( len(usno_line) > 0 ):
                    if re_comm.match(usno_line):
                        continue
                    elif re_html.search(usno_line):
                        error_msg = "Problem querying USNO database.  Returned an HTML FILE: \n %s" %usno_out
                        log.error(error_msg)
                        raise Exception, error_msg
                    
                    usno_star_splitted = usno_line.split()
                    
                    # Error control
                    if len(usno_star_splitted) != 11:
                        error_msg = "Unexpected USNO query response!  Length %d" %len(usno_star_splitted)
                        log.error(error_msg)
                        log.error("Response was:")
                        if( len(usno_catalogue) < 50 ):
                            for line in usno_catalogue:
                                log.error(line)
                        else:
                            log.error("too long to print.  First non-comment line:")
                            log.error(usno_line)
                        raise Exception, error_msg
                        
                    fields = {}
                    
                    # mapping
                    id = usno_star_splitted[10]
                    id = (int)(re.sub("\D", "", id))
                    fields['id'] = id
                    fields['ra'] = float(usno_star_splitted[0])
                    fields['dec'] = float(usno_star_splitted[1])
                    fields['ra_err'] = float(usno_star_splitted[4])/3600000.0
                    fields['dec_err'] = float(usno_star_splitted[5])/3600000.0
                    
                    if( usno_star_splitted[6] == "null" ):
                        if( usno_star_splitted[7] != "null" ):
                            fields['mag_b'] = float(usno_star_splitted[7])
                            fields['mag_err_b'] = 99.0
                        else:
                            fields['mag_b'] = 99.0
                            fields['mag_err_b'] = 99.0
                    else:
                        if( usno_star_splitted[7] == "null" ):
                            fields['mag_b'] = float(usno_star_splitted[6])
                            fields['mag_err_b'] = 99.0
                        else:
                            fields['mag_b'] = ( float(usno_star_splitted[6]) + float(usno_star_splitted[7]) ) / 2.0
                            fields['mag_err_b'] = abs( float(usno_star_splitted[6]) - float(usno_star_splitted[7]) )
                    if( usno_star_splitted[8] == "null" ):
                        if( usno_star_splitted[9] != "null" ):
                            fields['mag_r'] = float(usno_star_splitted[9])
                            fields['mag_err_r'] = 99.0
                        else:
                            fields['mag_r'] = 99.0
                            fields['mag_err_r'] = 99.0
                    else:
                        if( usno_star_splitted[9] == "null" ):
                            fields['mag_r'] = float(usno_star_splitted[8])
                            fields['mag_err_r'] = 99.0
                        else:
                            fields['mag_r'] = ( float(usno_star_splitted[8]) + float(usno_star_splitted[9]) ) / 2.0
                            fields['mag_err_r'] = abs( float(usno_star_splitted[8]) - float(usno_star_splitted[9]) )
                    star_tmp = model.USNO(**fields)
                    
                    # Add object to the list
                    self.objects.append(star_tmp)
                    
    
    def usno_load(self, sky_region):
        log.debug("Loading USNO data from local DB...")
        self.objects = []
        self.query = model.session.query(model.USNO).filter(model.USNO.ra > sky_region['ra']['min'],
                                                            model.USNO.ra < sky_region['ra']['max'],
                                                            model.USNO.dec > sky_region['dec']['min'],
                                                            model.USNO.dec < sky_region['dec']['max'])
        self.objects = self.query.all()
        
        # DETACH HERE so that we don't inadvertantly update the USNO table (e.g. when removing extinction)
        for obj in self.objects:
            model.session.expunge(obj)
            
        log.debug( "Loaded %d USNO objects from the local DB" %len(self.objects) )
        
        
    # SDSS catalog methods
    
    # Sync local db with sdss stars db at sky_region
    def sdss_load_external(self, sky_region, config):
        
        # query code is from sdss.py in pixelsim
        log.debug("Query to SDSS external data base...")
        import urllib
#        us_url = 'http://skyserver.sdss3.org/dr8/en/tools/search/x_sql.asp' #http://cas.sdss.org/astro/en/tools/search/x_sql.asp'
#        uk_url = 'http://www.sdss.org.uk/dr7/en/tools/search/x_sql.asp'
        us_url =  config['urls']['sdss'][0]
        uk_url =  config['urls']['sdss'][1]
        primary_url = us_url
        alternative_url = us_url # don't use another one because there is only one for dr8.
        
        fmt='csv'
        re_pass = re.compile("null|NONLEGACY")
        re_nogood = re.compile("GALAXY|QSO|STAR_CARBON|STAR_BROWN_DWARF|STAR_CATY_VAR|STAR_WHITE_DWARF")

        sky_regions = []
        if sky_region['ra']['max'] > sky_region['ra']['min']:
            sky_regions.append(sky_region)
        else:
            sky_region1 = dict(ra=dict(min=sky_region['ra']['min'],max=360.0),dec=dict(min=sky_region['dec']['min'],max=sky_region['dec']['max']))
            sky_region2 = dict(ra=dict(min=0.0,max=sky_region['ra']['max']),dec=dict(min=sky_region['dec']['min'],max=sky_region['dec']['max']))
            sky_regions.append(sky_region1)
            sky_regions.append(sky_region2)
            
        for sky_reg in sky_regions:
        
            sqlcmd =("select P.ra,P.dec,P.psfMag_u,P.psfMagErr_u,P.psfMag_g,P.psfMagErr_g,P.psfMag_r,P.psfMagErr_r,P.psfMag_i,P.psfMagErr_i,P.psfMag_z,P.psfMagErr_z," + 
                     "S.ObjType,P.objID,P.probPSF_u,P.probPSF_g,P.probPSF_r,P.probPSF_i,P.probPSF_z," + 
                     "P.raErr,P.decErr,P.type,P.probPSF " + 
                     "from star as P left outer join SpecObj as S on P.objID = S.BestObjID " + 
                     "where P.ra between %f and %f and P.dec between %f and %f and P.clean=1" %(sky_reg['ra']['min'], sky_reg['ra']['max'], sky_reg['dec']['min'], sky_reg['dec']['max']))
            params = urllib.urlencode({'cmd': sqlcmd, 'format': fmt})
            
            # Call to primary mirror
            sdss_out = (urllib.urlopen(primary_url + '?%s' % params)).read()
            log.debug("...selection received.")
            
            sdss_catalogue = sdss_out.split('\n')
            
            # Check primary
            #log.debug("SDSS Primary mirror header response: %s" %sdss_catalogue[0].split(','))
            if len(sdss_catalogue[0].split(',')) != 23:
                log.warning("Primary SDSS mirror failed. Trying with alternative...")
                # Call to alternative
                sdss_out = (urllib.urlopen(alternative_url + '?%s' % params)).read()
                sdss_catalogue = sdss_out.split('\n')
                log.debug("SDSS Alternative mirror header response: %s" %sdss_catalogue[0].split(','))
                if len(sdss_catalogue[0].split(',')) != 23:
                    error_msg = "Unexpected SDSS query header for both primary and alternative mirrors!"
                    log.error(error_msg)
                    raise Exception, error_msg
            
            if len (sdss_catalogue) == 1:
                # No objects found
                error_msg = "No objects found in SDSS catalogue!"
                log.error(error_msg)
                raise Exception, error_msg
            else:
                # remove first and last line 
                sdss_catalogue.remove(sdss_catalogue[0])
                sdss_catalogue.remove(sdss_catalogue[-1])
                log.debug("Found %d stars in SDSS catalogue" %len(sdss_catalogue))
            
            # Main loop
            
            for sdss_star in sdss_catalogue:
                sdss_star_splitted = sdss_star.split(',')
            
                # Error control
                if len(sdss_star_splitted) != 23:
                    error_msg = "Unexpected SDSS query element response!"
                    log.error(error_msg)
                    log.error("Response was:")
                    if( len(sdss_catalogue) < 50 ):
                        for line in sdss_catalogue:
                            log.error(line)
                    else:
                        log.error("too long to print.  First line:")
                        log.error(sdss_star)
                    raise Exception, error_msg
            
                # throw away any objects that are spectrally determined 
                # to be not what we want: qsos, white dwarfs, carbon stars, 
                # brown dwarfs, galaxies.
                specType = sdss_star_splitted[12]
                good = True
                #if re.search("null", specType) is not None or re.search("NONLEGACY", specType) is not None:
                if re_pass.search(specType):
                    pass
                #elif re.search("GALAXY", specType) is not None or re.search("QSO", specType) is not None or re.search("STAR_CARBON", specType) is not None or re.search("STAR_BROWN_DWARF", specType) is not None or re.search("STAR_CATY_VAR", specType) is not None or re.search("STAR_WHITE_DWARF", specType) is not None:
                elif re_nogood.search(specType):
                    good = False
                elif isinstance( specType, int ):
                    typeNum = int(specType)
                    if typeNum == 0 or typeNum == 1 or typeNum == 3 or typeNum == 4 or typeNum == 6 or typeNum == 8:
                        good = False
                if not good:
                    continue
                
                # Be careful! Sdss magnitudes do not match DES nor PAU!
                # Use for compute star types.
                
                fields = {}
                # mapping
                
                fields['id'] = int(sdss_star_splitted[13]) 
                fields['ra'] = float(sdss_star_splitted[0])
                fields['dec'] = float(sdss_star_splitted[1])
                fields['dec_err'] = 0.0
                fields['mag_u'] = float(sdss_star_splitted[2]) # u mag
                fields['mag_err_u'] = float(sdss_star_splitted[3])
                fields['mag_g'] = float(sdss_star_splitted[4]) # g mag
                fields['mag_err_g'] = float(sdss_star_splitted[5])
                fields['mag_r'] = float(sdss_star_splitted[6]) # r mag
                fields['mag_err_r'] = float(sdss_star_splitted[7])
                fields['mag_i'] = float(sdss_star_splitted[8]) # i mag
                fields['mag_err_i'] = float(sdss_star_splitted[9])
                fields['mag_z'] = float(sdss_star_splitted[10]) # z mag            
                fields['mag_err_z'] = float(sdss_star_splitted[11])
                
                fields['u_s'] = float(sdss_star_splitted[14])
                fields['g_s'] = float(sdss_star_splitted[15])
                fields['r_s'] = float(sdss_star_splitted[16])
                fields['i_s'] = float(sdss_star_splitted[17])
                fields['z_s'] = float(sdss_star_splitted[18])
                
                # N/A
                fields['ra_err'] = float(sdss_star_splitted[19])
                fields['dec_err'] = float(sdss_star_splitted[20])
                fields['type'] = int(sdss_star_splitted[21])
                fields['prob_psf'] = float(sdss_star_splitted[22])
                
                star_tmp = model.SDSS(**fields)
                
                # Add object to the list
                self.objects.append(star_tmp)
                
    
    def sdss_load(self, sky_region):
        log.debug("Loading SDSS data from local DB...")
        self.objects = []
        self.query = model.session.query(model.SDSS).filter(model.SDSS.ra > sky_region['ra']['min'],
                                                            model.SDSS.ra < sky_region['ra']['max'],
                                                            model.SDSS.dec > sky_region['dec']['min'],
                                                            model.SDSS.dec < sky_region['dec']['max'])
        self.objects = self.query.all()
        
        # DETACH HERE so that we don't inadvertantly update the SDSS table (e.g. when removing extinction)
        for obj in self.objects:
            model.session.expunge(obj)
        
        log.debug( "Loaded %d SDSS objects from the local DB" %len(self.objects) )
        
    
    def removeExtinction(self, instrument):
        
        log.debug( "Removing extinction!" )
        
        # load schlegel maps
        hdulist_n = pyfits.open(resource_filename('paudm.resources','extinct/SFD_dust_4096_ngp.fits'))
        SFD_north = hdulist_n[0].data
        hdulist_s = pyfits.open(resource_filename('paudm.resources','extinct/SFD_dust_4096_sgp.fits'))
        SFD_south = hdulist_s[0].data
        
        # assuming sdss filter constants from Schlafly & Finkbeiner 2010
        # and extrapolating to the narrow bands.
        # the outermost values are linear interpolations, to prevent bounds errors.
        sdss_bb_wvls = [122.7, 354.3, 477.0, 623.1, 762.5, 913.4, 1064.3]
        sdss_extinct = [5.175, 4.239, 3.303, 2.285, 1.698, 1.263, 0.828]
        extinction_coeffs_func = interp1d(sdss_bb_wvls, sdss_extinct, kind='cubic')
        
        
        for object in self.objects:
            
            # find galactic longitude and latitude
            # algorith from 
            # http://www.astro.rug.nl/software/kapteyn/celestialbackground.html
            ra_rad = object.ra*math.pi/180.
            dec_rad = object.dec*math.pi/180.
            b = math.asin(-0.867666135683*math.cos(dec_rad)*math.cos(ra_rad) - 0.198076389613*math.cos(dec_rad)*math.sin(ra_rad) + 0.455983794521*math.sin(dec_rad))
            cosl = (-0.054875539396*math.cos(dec_rad)*math.cos(ra_rad) - 0.873437104728*math.cos(dec_rad)*math.sin(ra_rad) - 0.48383499177*math.sin(dec_rad))/math.cos(b)
            if( cosl > 1.0 and cosl < 1.0000000001 ):
                cosl = 1.0;
            if( math.fabs(cosl) > 1. ):
                error_msg = "starMaker.removeExtinction: conversion to galactic coordinates failed!"
                log.error(error_msg)
                raise Exception, error_msg
            
            sinl = (0.494109453628*math.cos(dec_rad)*math.cos(ra_rad) - 0.444829594298*math.cos(dec_rad)*math.sin(ra_rad) + 0.7469822487*math.sin(dec_rad))/math.cos(b)
            if( sinl > 1.0 and sinl < 1.0000000001 ):
                sinl = 1.0;
            if( math.fabs(sinl) > 1. ):
                error_msg = "starMaker.removeExtinction: conversion to galactic coordinates failed!"
                log.error(error_msg)
                raise Exception, error_msg
            
            l = math.atan2(sinl,cosl);
            if( l < 0. ):
                l = l + 2.*math.pi;
            
            # find x and y on the maps
            if b >= 0.:
                n = 1
            else:
                n = -1
            x = 2048. * math.sqrt(1-n*math.sin(b))*math.cos(l) + 2047.5
            y = -2048. * n * math.sqrt(1-n*math.sin(b))*math.sin(l) + 2047.5
            
            # look up the value of that pixel in the dust map
            x_int = math.floor(x+0.5)
            y_int = math.floor(y+0.5)
            if n == 1:
                extinction = SFD_north[y_int-1, x_int-1];
            else:
                extinction = SFD_south[y_int-1, x_int-1];
                
                
            if( type(object) == model.SDSS ):
                # and remove the extinction 
                # assuming sdss filter constants from Schlafly & Finkbeiner 2010
                object.mag_u -= 4.239*extinction
                object.mag_g -= 3.303*extinction
                object.mag_r -= 2.285*extinction
                object.mag_i -= 1.698*extinction
                object.mag_z -= 1.263*extinction
            elif( type(object) == model.Truth_Object ):
                for filter_name in instrument['BROAD_FILTERS'] + instrument['NARROW_FILTERS']:
                    extinction_coeff = extinction_coeffs_func(instrument['FILTER_ASSOCIATIONS'][filter_name]['WAVELEN'])
                    magname = "mag_" + "%s" %filter_name
                    newmag = object.__getattribute__(magname)
                    newmag -= extinction_coeff*extinction
                    setattr(object, magname, newmag)
            elif( type(object) == model.USNO ):
                # just use g and r values, its close enough, since the USNO photometry sucks.
                object.mag_b -= 3.303*extinction
                object.mag_r -= 2.285*extinction
            else:
                error_msg = "removeExtinction: Trying to run on a non-SDSS or Truth object!  Fix me!"
                log.error(error_msg)
                raise Exception, error_msg
            
            
    
    # extinct the ideal objects to make realistic observations
    def addExtinction(self, instrument, detection_filter=None):
        
        log.debug( "Adding extinction!" )
        
        # load schlegel maps
        hdulist_n = pyfits.open(resource_filename('paudm.resources',"extinct/SFD_dust_4096_ngp.fits"))
        SFD_north = hdulist_n[0].data
        hdulist_s = pyfits.open(resource_filename('paudm.resources',"extinct/SFD_dust_4096_sgp.fits"))
        SFD_south = hdulist_s[0].data
        
        # assuming sdss filter constants from Schlafly & Finkbeiner 2010
        # and extrapolating to the narrow bands.
        # the outermost values are linear interpolations, to prevent bounds errors.
        sdss_bb_wvls = [122.7, 354.3, 477.0, 623.1, 762.5, 913.4, 1064.3]
        sdss_extinct = [5.175, 4.239, 3.303, 2.285, 1.698, 1.263, 0.828]
        extinction_coeffs_func = interp1d(sdss_bb_wvls, sdss_extinct, kind='cubic')
        
        for object in self.objects:
            
            # find galactic longitude and latitude
            # algorith from 
            # http://www.astro.rug.nl/software/kapteyn/celestialbackground.html
            ra_rad = object.ra*math.pi/180.
            dec_rad = object.dec*math.pi/180.
            b = math.asin(-0.867666135683*math.cos(dec_rad)*math.cos(ra_rad) - 0.198076389613*math.cos(dec_rad)*math.sin(ra_rad) + 0.455983794521*math.sin(dec_rad))
            cosl = (-0.054875539396*math.cos(dec_rad)*math.cos(ra_rad) - 0.873437104728*math.cos(dec_rad)*math.sin(ra_rad) - 0.48383499177*math.sin(dec_rad))/math.cos(b)
            if( cosl > 1.0 and cosl < 1.0000000001 ):
                cosl = 1.0;
            if( math.fabs(cosl) > 1. ):
                error_msg = "starMaker.removeExtinction: conversion to galactic coordinates failed!"
                log.error(error_msg)
                raise Exception, error_msg
                
            sinl = (0.494109453628*math.cos(dec_rad)*math.cos(ra_rad) - 0.444829594298*math.cos(dec_rad)*math.sin(ra_rad) + 0.7469822487*math.sin(dec_rad))/math.cos(b)
            if( sinl > 1.0 and sinl < 1.0000000001 ):
                sinl = 1.0;
            if( math.fabs(sinl) > 1. ):
                error_msg = "starMaker.removeExtinction: conversion to galactic coordinates failed!"
                log.error(error_msg)
                raise Exception, error_msg
                
            l = math.atan2(sinl,cosl)
            if( l < 0. ):
                l = l + 2.*math.pi;
                
            # find x and y on the maps
            if b >= 0.:
                n = 1
            else:
                n = -1
            x = 2048. * math.sqrt(1-n*math.sin(b))*math.cos(l) + 2047.5
            y = -2048. * n * math.sqrt(1-n*math.sin(b))*math.sin(l) + 2047.5
            
            # look up the value of that pixel in the dust map
            x_int = math.floor(x+0.5)
            y_int = math.floor(y+0.5)
            if n == 1:
                extinction = SFD_north[y_int-1, x_int-1];
            else:
                extinction = SFD_south[y_int-1, x_int-1];
                
            
            # and add the extinction 
            
            # if the object is a Truth_Object or an SDSS DB object
            if type(object) == model.Truth_Object:
                for filter_name in instrument['BROAD_FILTERS'] + instrument['NARROW_FILTERS']:
                    extinction_coeff = extinction_coeffs_func(instrument['FILTER_ASSOCIATIONS'][filter_name]['WAVELEN'])
                    magname = "mag_" + "%s" %filter_name
                    newmag = object.__getattribute__(magname)
                    newmag += extinction_coeff*extinction
                    setattr(object, magname, newmag)
            
            elif type(object) == model.SDSS:
                object.mag_u += 4.239*extinction
                object.mag_g += 3.303*extinction
                object.mag_r += 2.285*extinction
                object.mag_i += 1.698*extinction
                object.mag_z += 1.263*extinction
                
            elif type(object) == model.USNO:
                # just use g and r values, its close enough, since the USNO photometry sucks.
                object.mag_b += 3.303*extinction
                object.mag_r += 2.285*extinction
            
            # if the object is a Detection
            elif type(object) == model.Detection:
                if detection_filter is None:
                    log.warning( "No detection filter specified in addExtinction!" )
                    break
                extinction_coeff = extinction_coeffs_func(instrument['FILTER_ASSOCIATIONS'][detection_filter]['WAVELEN'])
                fluxname = "flux_auto"
                newflux = object.__getattribute__(fluxname)
                newflux /= 10.0**(0.4*extinction_coeff*extinction)
                setattr(object, fluxname, newflux)
                
            else:
                log.warning( "Object type %s not handled in addExtinction!" )
                break
                
    
    # returns an sdss-like error appropriate for a given mag and filter
    # the data are taken from predictions in Gunn's SDSS camera paper:
    # http://iopscience.iop.org/1538-3881/116/6/3040/fulltext/
    # Table 5, plus an ubercal-like systematic error.
    def sdss_error_from_mag( self, mag, filter ):
        
        # These systematic errors are meant to be ubercal-appropriate
        systematics = dict()
        systematics['u'] = 0.02
        systematics['g'] = 0.01
        systematics['r'] = 0.01
        systematics['i'] = 0.01
        systematics['z'] = 0.01
        
        # This info is from the Gunn paper cited above this function
        mags = [17.0, 17.5, 18.0, 18.5, 19.0, 19.5, 20.0, 20.5, 21.0, 21.5, 22.0, 22.5, 23.0, 23.5, 24.0, 24.5, 25.0]
        errors = dict()
        errors['u'] = [0.00623441, 0.00803213, 0.0104712, 0.0138696, 0.0187617, 0.025974, 0.0369004, 0.0537634, 0.0806452, 0.121951, 0.188679, 0.294118, 0.454545, 0.714286, 1.11111, 1.66667, 2.5]
        errors['g'] = [0.00233536, 0.00299043, 0.00386399, 0.00505306, 0.00672495, 0.00914077, 0.0127389, 0.0182482, 0.026738, 0.0401606, 0.0609756, 0.0943396, 0.147059, 0.232558, 0.357143, 0.588235, 0.909091]
        errors['r'] = [0.00237982, 0.00307692, 0.00403063, 0.00537346, 0.00732064, 0.0102249, 0.0146628, 0.0215517, 0.0323625, 0.049505, 0.0763359, 0.119048, 0.185185, 0.294118, 0.454545, 0.714286, 1.11111]
        errors['i'] = [0.00295683, 0.00392157, 0.00531067, 0.0073692, 0.0105042, 0.0153374, 0.0229358, 0.0348432, 0.0537634, 0.0833333, 0.131579, 0.204082, 0.322581, 0.526316, 0.833333, 1.25, 2.0]
        errors['z'] = [0.00805153, 0.011534, 0.0169492, 0.0254453, 0.0389105, 0.0598802, 0.0934579, 0.147059, 0.232558, 0.37037, 0.588235, 0.909091, 1.42857, 2.5, 3.33333, 5.0, 10.0]
        
        if filter not in errors:
            # we're trying to do a narrow band, probably... or Y!
            # which broad band is the narrow band closest to?  this will give an error that does not include sed-matching problems.
            wvl = 0
            try:
                wvl = instrument['FILTER_ASSOCIATIONS'][filter]['WAVELEN']
            except KeyError:
                error_msg = "Error, trying to assign an sdss-like error when the filter is not known: %s" %filter
                log.error(error_msg)
                raise Exception, error_msg
                
            best_dist = 999999.
            close_broadband = ''
            
            # Loop over the broad band filters with error information (All Broad bands except Y)
            for bfilter in errors.keys():
                broad_wvl = instrument['FILTER_ASSOCIATIONS'][bfilter]['WAVELEN']
                wvl_dist = abs(broad_wvl - wvl)
                if( wvl_dist < best_dist ):
                    best_dist = wvl_dist
                    close_broadband = bfilter
            filter = close_broadband
            
        statistical_error_func = interp1d( mags, errors[filter] )
        
        if mag < 17.0:
            mag = 17.0
        elif mag > 25.0:
            mag = 25.0
            
        # the error function really really returns N/S, so flux_error/flux rather than mag error, hence the 1.0857
        mag_stat_error = 1.0857*statistical_error_func(mag)
        error = math.sqrt(mag_stat_error**2 + systematics[filter]**2)
        
        return error
        
    
    def sdss_mag_from_error( self, error, filter ):
        
        # These systematic errors are meant to be ubercal-appropriate
        systematics = dict()
        systematics['u'] = 0.02
        systematics['g'] = 0.01
        systematics['r'] = 0.01
        systematics['i'] = 0.01
        systematics['z'] = 0.01
        
        # This info is from the Gunn paper cited above sdss_error_from_mag
        mags = [17.0, 17.5, 18.0, 18.5, 19.0, 19.5, 20.0, 20.5, 21.0, 21.5, 22.0, 22.5, 23.0, 23.5, 24.0, 24.5, 25.0]
        errors = dict()
        errors['u'] = [0.00623441, 0.00803213, 0.0104712, 0.0138696, 0.0187617, 0.025974, 0.0369004, 0.0537634, 0.0806452, 0.121951, 0.188679, 0.294118, 0.454545, 0.714286, 1.11111, 1.66667, 2.5]
        errors['g'] = [0.00233536, 0.00299043, 0.00386399, 0.00505306, 0.00672495, 0.00914077, 0.0127389, 0.0182482, 0.026738, 0.0401606, 0.0609756, 0.0943396, 0.147059, 0.232558, 0.357143, 0.588235, 0.909091]
        errors['r'] = [0.00237982, 0.00307692, 0.00403063, 0.00537346, 0.00732064, 0.0102249, 0.0146628, 0.0215517, 0.0323625, 0.049505, 0.0763359, 0.119048, 0.185185, 0.294118, 0.454545, 0.714286, 1.11111]
        errors['i'] = [0.00295683, 0.00392157, 0.00531067, 0.0073692, 0.0105042, 0.0153374, 0.0229358, 0.0348432, 0.0537634, 0.0833333, 0.131579, 0.204082, 0.322581, 0.526316, 0.833333, 1.25, 2.0]
        errors['z'] = [0.00805153, 0.011534, 0.0169492, 0.0254453, 0.0389105, 0.0598802, 0.0934579, 0.147059, 0.232558, 0.37037, 0.588235, 0.909091, 1.42857, 2.5, 3.33333, 5.0, 10.0]
        
        if filter not in errors:
            # we're trying to do a narrow band, probably...
            # which broad band is the narrow band closest to?  this will give an error that does not include sed-matching problems.
            wvl = 0
            try:
                wvl = instrument['FILTER_ASSOCIATIONS'][filter]['WAVELEN']
            except KeyError:
                error_msg = "Error, trying to assign an sdss-like mag when the filter is not known: %s" %filter
                log.error(error_msg)
                raise Exception, error_msg
                
            best_dist = 999999.
            close_broadband = ''
            for bfilter in instrument['BROAD_FILTERS']:
                broad_wvl = instrument['FILTER_ASSOCIATIONS'][bfilter]['WAVELEN']
                wvl_dist = abs(broad_wvl - wvl)
                if( wvl_dist < best_dist ):
                    best_dist = wvl_dist
                    close_broadband = bfilter
            filter = close_broadband
            
        # now add the statistical and systematic errors
        for i in range(len(errors[filter])):
            errors[filter][i] = math.sqrt(errors[filter][i]**2 + systematics[filter]**2)
            
        mag = 99.9
        if( error < errors[filter][0] ):
            mag = mags[0]
        elif( error > errors[filter][-1] ):
            mag = mags[-1]
        else:
            mag_func = interp1d( errors[filter], mags )
            mag = float(mag_func(error))
            
        return mag
        
    
    def truth_to_sdss(self, sdss_filter_system, pau_filter_system, star_seds):
        
        # remove the galactic extinction
        self.removeExtinction()
        
        # translate the SEDs from floats into strings
        sed_string_indexing = star_seds.seds.keys()
        
        index=0
        sdss_objects = []
        for s in range( len(self.objects) ):
            
            fields = {}
            fields['id'] = self.objects[s].id
            fields['ra'] = self.objects[s].ra
            fields['dec'] = self.objects[s].dec
            fields['ra_err'] = 0.000028 # 0.1 arcsec.  not used anywhere so far...
            fields['dec_err'] = 0.000028
            fields['type'] = 6
            fields['prob_psf'] = 1
            fields['u_s'] = 1
            fields['g_s'] = 1
            fields['r_s'] = 1
            fields['i_s'] = 1
            fields['z_s'] = 1
            star_tmp = model.SDSS(**fields)
            
            norm_mag = self.objects[s].mag_i
            for filter in( 'u', 'g', 'r', 'i', 'z' ):
                # convolve with the sed_type to get the SDSS magnitudes
                flux = sdss_filter_system.filters[filter].sed_to_flux(star_seds.seds[sed_string_indexing[(int)(self.objects[s].sed_type)]], pau_filter_system.filters['i'], norm_mag, True)
                newmag = constants.mag_zp - 2.5*math.log10(flux)
                magname = "mag_" + filter
                setattr(star_tmp, magname, newmag)
                magerr = self.sdss_error_from_mag( newmag, filter )
                magname = "mag_err_" + filter
                setattr(star_tmp, magname, magerr)
                
            sdss_objects.append(star_tmp)
            
            index = index+1
            
        self.objects = sdss_objects
        
        # add extinction to the magnitudes
        self.addExtinction()
    
    def sdss_to_detection(self, filter, sdss_filter_system, pau_filter_system, star_seds):
        
        # remove the galactic extinction
        self.removeExtinction()
        
        # estimate which stellar types the stars are
        best_seds = self.sdss_assign_seds(sdss_filter_system, star_seds)
        
        if( len(best_seds) != len(self.objects) ):
            error_msg = "best sed list is not the proper length: %d != %d !" %(len(best_seds), len(self.objects))
            log.error(error_msg)
            
        index=0
        detection_objects = []
        for s in range( len(self.objects) ):
            
            # mapping from sdss object to detection object
            # (which has only one magnitude:  that of the chip in question)
            fields = {}
            fields['id'] = self.objects[s].id
            fields['ra'] = self.objects[s].ra
            fields['dec'] = self.objects[s].dec
            
            # this is not actually correct for the faintest objects where pogson mags are different from asinh mags...
            # but we should not be using those in our analyses anyway!
            fields['flux_auto'] = pau_filter_system.filters[filter].sed_to_flux(star_seds.seds[best_seds[s]], sdss_filter_system.filters['i'], self.objects[s].mag_i, True)
            fields['flux_err_auto'] = math.fabs(10**(self.objects[s].mag_err_i/2.5)*fields['flux_auto'] - fields['flux_auto'])
            
            # fill in all flux fields with the same info
            fields['flux_psf'] = fields['flux_auto']
            fields['flux_err_psf'] = fields['flux_err_auto']
            fields['flux_model'] = fields['flux_auto']
            fields['flux_err_model'] = fields['flux_err_auto']
            
            star_tmp = model.Detection(**fields)
            detection_objects.append(star_tmp)
            
            index = index+1
            
        self.objects = detection_objects
        
        # add extinction to the magnitudes
        self.addExtinction( filter )
        
    def truth_to_detection(self, filter):
        
        detection_objects = []
        for s in range( len(self.objects) ):
            
            # mapping from sdss object to detection object
            # (which has only one magnitude:  that of the chip in question)
            fields = {}
            fields['id'] = self.objects[s].id
            fields['ra'] = self.objects[s].ra
            fields['dec'] = self.objects[s].dec
            
            # this is not actually correct for the faintest objects where pogson mags are different from asinh mags...
            # but we should not be using those in our analyses anyway!
            mag_name = "mag_" + "%s" %filter
            mag = self.objects[s].__getattribute__(mag_name)
            fields['flux_auto'] = 10.0**((constants.mag_zp-mag)/2.5)
            mag_err = self.sdss_error_from_mag( mag, filter )
            fields['flux_err_auto'] = math.fabs(10.0**(mag_err/2.5)*fields['flux_auto'] - fields['flux_auto'])
            
            # fill in all flux fields with the same info
            fields['flux_psf'] = fields['flux_auto']
            fields['flux_err_psf'] = fields['flux_err_auto']
            fields['flux_model'] = fields['flux_auto']
            fields['flux_err_model'] = fields['flux_err_auto']
            
            star_tmp = model.Detection(**fields)
            detection_objects.append(star_tmp)
            
        self.objects = detection_objects
    
    def sdss_assign_seds(self, filter_system, star_seds):
        
        colors = {}
        for sed in star_seds.seds.keys():
            colors[sed] = {}
        
        # first, for every sed calculate the sdss colors.
        for sed in star_seds.seds.keys():
            uflux = filter_system.filters['u'].sed_to_flux( star_seds.seds[sed], filter_system.filters['i'], 19.0, True )
            gflux = filter_system.filters['g'].sed_to_flux( star_seds.seds[sed], filter_system.filters['i'], 19.0, True )
            rflux = filter_system.filters['r'].sed_to_flux( star_seds.seds[sed], filter_system.filters['i'], 19.0, True )
            iflux = filter_system.filters['i'].sed_to_flux( star_seds.seds[sed], filter_system.filters['i'], 19.0, True )
            zflux = filter_system.filters['z'].sed_to_flux( star_seds.seds[sed], filter_system.filters['i'], 19.0, True )
            
            colors[sed]['u_minus_g'] = 2.5*math.log10(gflux/uflux)
            colors[sed]['g_minus_r'] = 2.5*math.log10(rflux/gflux)
            colors[sed]['r_minus_i'] = 2.5*math.log10(iflux/rflux)
            colors[sed]['i_minus_z'] = 2.5*math.log10(zflux/iflux)
            
        # now, for every star, see which seds have consistent colors.
        no_match3 = 0
        no_match5 = 0
        
        best_seds = []
        for star in self.objects:
            
            # add 1% systematic error to each mag (ubercal error estimate)
            g_minus_r = star.mag_g - star.mag_r
            g_minus_r_error = math.sqrt(star.mag_err_g*star.mag_err_g + star.mag_err_r*star.mag_err_r + 0.0002)
            r_minus_i = star.mag_r - star.mag_i
            r_minus_i_error = math.sqrt(star.mag_err_r*star.mag_err_r + star.mag_err_i*star.mag_err_i + 0.0002)
            i_minus_z = star.mag_i - star.mag_z
            i_minus_z_error = math.sqrt(star.mag_err_i*star.mag_err_i + star.mag_err_z*star.mag_err_z + 0.0002)
            # add a huge term to the u error to account 
            # for the red leak in the u band
            # this is not an exact correction, since the effect 
            # depends on the instrumental r-i.  see the sdss website.
            u_minus_g = star.mag_u - star.mag_g
            u_minus_g_error = math.sqrt(star.mag_err_u*star.mag_err_u + star.mag_err_g*star.mag_err_g + 0.09 )
            
            possibleSEDs = []
            nsigma = 0
            worstcolor = 0
            worstNsigma = 1e9
            bestSED = None
            for sed in star_seds.seds.keys():
                # don't use u-g?
                ug_sigma = 0
                #ug_sigma = math.fabs(colors[sed]['u_minus_g'] - u_minus_g)/u_minus_g_error
                gr_sigma = math.fabs(colors[sed]['g_minus_r'] - g_minus_r)/g_minus_r_error
                ri_sigma = math.fabs(colors[sed]['r_minus_i'] - r_minus_i)/r_minus_i_error
                # TEST: don't use z???
                iz_sigma = math.fabs(colors[sed]['i_minus_z'] - i_minus_z)/i_minus_z_error
                #iz_sigma = 0
                if max(ug_sigma,gr_sigma,ri_sigma,iz_sigma) < math.fabs(worstNsigma):
                    worstNsigma = max(ug_sigma,gr_sigma,ri_sigma,iz_sigma)
                    bestSED = sed
                    if worstNsigma == ug_sigma:
                        worstcolor = 1
                        worstNsigma = (colors[sed]['u_minus_g'] - u_minus_g)/u_minus_g_error
                    elif worstNsigma == gr_sigma:
                        worstcolor = 2
                        worstNsigma = (colors[sed]['g_minus_r'] - g_minus_r)/g_minus_r_error
                    elif worstNsigma == ri_sigma:
                        worstcolor = 3
                        worstNsigma = (colors[sed]['r_minus_i'] - r_minus_i)/r_minus_i_error
                    elif worstNsigma == iz_sigma:
                        worstcolor = 4
                        worstNsigma = (colors[sed]['i_minus_z'] - i_minus_z)/i_minus_z_error 
                        
            best_seds.append(bestSED)
            
            if math.fabs(worstNsigma) > 5.:
                no_match5 = no_match5+1
            if math.fabs(worstNsigma) > 3.:
                no_match3 = no_match3+1
                
        log.debug("%d (%d) with no SED matches within 3 (5) sigma" %(no_match3, no_match5))
        
        return best_seds
    
    
    # GENERAL Catalogue methods
    def infer_query_area(self):
        """
        Search over the query what are the field limits
        """
        
        field_range = {'ra' : {'min': 360, 'max': 0},
                       'dec': {'min': 90, 'max': -90}}
        
        
        for element in self.query:
            if element.ra < field_range['ra']['min']:
                field_range['ra']['min'] = element.ra
            if element.ra > field_range['ra']['max']:
                field_range['ra']['max'] = element.ra
            if element.dec < field_range['dec']['min']:
                field_range['dec']['min'] = element.dec
            if element.dec > field_range['dec']['max']:
                field_range['dec']['max'] = element.dec
            
        return field_range
    
    
    def dump2ascii(self, columns, filename, input_type='query', separator='\t', ):
        """
        Export the catalogue to an ascii file
        """
        fd = open(filename, 'w')
        
        log.info("Exporting catalogue from " + input_type + "...")
        
        # Header
        fd.write('#')
        for parameter in columns:
            fd.write(parameter + separator)
        fd.write('\n')
        
        # Looping queries
        if input_type == 'query':
            for element in self.query:
                for parameter in columns:
                    fd.write(str(element.__getattribute__(parameter)) + separator)
                fd.write('\n')
        
        # Looping objects
        if input_type == 'objects':
            for element in self.objects:
                for parameter in columns:
                    fd.write(str(element.__getattribute__(parameter)) + separator)
                fd.write('\n')
        
        fd.close()
        log.info("Catalogue correctly exported to file: " + filename)
        
    
    def dump2ROOT(self, columns, filename, separator='\t', ):
        """
        Export the catalogue to a ROOT file
        """
        pass # UNUSED! (Erase if code is OBSOLETE)
        return -1
        """
        from ROOT import TFile,TTree
        from array import array
        
        fd = TFile(filename,'RECREATE')
        td = TTree('TPAU','Tree with catalog data')
        
        log.info("Creating a ROOT file from catalog")
        
        # Create branches
        n = len(columns)
        var = array('f', n*[0.] )
        for element in self.query:
            i = 0
            for parameter in columns:
                td.Branch(parameter,var[i],parameter+'/F')
                i = i + 1
        
        # Looping queries
        for element in self.query:
            i = 0
            for parameter in columns:
                var[i] = str(element.__getattribute__(parameter))
                i = i + 1
            td.Fill()
            
        fd.Write()
        fd.Close()
        log.info("Catalogue correctly exported to ROOT file: " + filename)
        """
    
    
    # GENERAL Plotting methods
    def plot_x(self, parameter, use_query=True, xlabel = "X Label", ylabel = "Y Label", title = "Title", filename = None ,format=None):
        """Plot a parameter given from the available catalogue.
        If range is not given, it is adjusted automatically.
        if filename is given, it is saved and stored. If not, it is displayed in the screen.
        If use_query is True, catalogue elements are loaded directly from the data base (self.query). Use if possible.
        If use_query is false, catalogue elements are loaded from memory (self.objects)
        """
        import matplotlib.pyplot as plt
        
        # Load parameter array
        x_array = []
        
        # Select input mode
        if use_query:
            iterator = self.query
        else:
            iterator = self.objects
        
        # Loop
        if format == None: 
            for db_object in iterator:
                x_array.append(db_object.__getattribute__(parameter))
        else:
            x_array = parameter
        
        # Create plot
        plt.plot(x_array)
        
        # Configure plot
        if xlabel == None:
            xlabel = ''
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.title(title)
        
        # Show or save
        if filename == None:
            plt.show()
        else:
            plt.savefig(filename)
        
    
    def plot_xy(self, parameter_x, parameter_y, type = '.', equal_axis = False, use_query=True, xlabel = None, ylabel = None, title = "Title", filename = None, color = 'b', range_x=None, range_y=None ,option=None, format=None):
        """Plot two parameters given from the available catalogue 
        If range is not given, it is adjusted automatically
        if filename is given, it is saved and stored. If not, it is displayed in the screen.
        Type can be dot or line.
        If use_query is True, catalogue elements are loaded directly from the data base (self.query). Use if possible.
        If use_query is false, catalogue elements are loaded from memory (self.objects)
        """
        import matplotlib.pyplot as plt
        
        #fig = plt.figure()
        plt.clf()
        
        # Load parameter array
        x_array = []
        y_array = []
        
        # Select input mode
        if use_query:
            iterator = self.query
        else:
            iterator = self.objects
        
        # Loop
        count = 0
        log.debug('Starting loop')
        if format == None:
            log.debug('Format is none')
            for db_object in iterator:
                if count % 100 == 0:
                    log.debug('Loading object ')
                valuex = db_object.__getattribute__(parameter_x)
                valuey = db_object.__getattribute__(parameter_y)
                
                if parameter_x == 'mag_auto':
                    image = model.session.query(model.Image).join(model.Detection).filter(model.Detection.id == db_object.id).one()
                    mosaic = model.session.query(model.Mosaic).join(model.Image).filter(model.Image.id == image.id).one()
                    zp_nightly = image.zp_nightly
                    ##if type(zp_nightly).__name__=='NoneType':
                        ##continue
                    airmass = mosaic.air_mass
                    valuex = valuex + zp_nightly*airmass
                if parameter_y == 'mag_auto':
                    image = model.session.query(model.Image).join(model.Detection).filter(model.Detection.id == db_object.id).one()
                    mosaic = model.session.query(model.Mosaic).join(model.Image).filter(model.Image.id == image.id).one()
                    zp_nightly = image.zp_nightly
                    ##if type(zp_nightly).__name__=='NoneType':
                        ##continue
                    airmass = mosaic.air_mass
                    valuey = valuey + zp_nightly*airmass
                
                x_array.append(valuex)
                y_array.append(valuey)
                count = count + 1
        else:
            x_array = parameter_x
            y_array = parameter_y
        
        # Create XY plot
        if option == 'heat':
            heatmap, xedges, yedges = np.histogram2d(x_array, y_array, bins=150)
            extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
            plt.imshow(heatmap, extent=extent)
        else:
            plt.plot(x_array, y_array, type, color = color, markersize=1)
        
        
        # Configure plot
        if xlabel == None:
            xlabel = ''
        plt.xlabel(xlabel)
        if ylabel == None:
            ylabel = ''
        plt.ylabel(ylabel)
        plt.title(title)
        if range_x != None and range_y != None:
            range_x = np.array(range_x)
            range_y = np.array(range_y)
            plt.axis(np.append(range_x,range_y))
        if equal_axis == True:
            plt.axis('equal')
            
        # Show or save
        if filename == None:
            plt.show()
        else:
            #plt.savefig(filename)
            plt.savefig(filename)
            
        plt.hold(False) 
    
    def plot_histogram(self, parameter, bins=10, use_query=True, xlabel=None, ylabel="Counts", range=None, title="Title", filename=None ,format = None):
        """Plot histogram of the parameter given from the available catalogue.
        If range is not given, it is adjusted automatically.
        If bins is not given, it is adjusted to 10.
        if filename is given, it is saved and stored. If not, it is displayed in the screen.
        If use_query is True, catalogue elements are loaded directly from the data base (self.query). Use if possible.
        If use_query is false, catalogue elements are loaded from memory (self.objects)
        """
        import matplotlib.pyplot as plt
        
        plt.figure()
        
        # Load parameter array
        #x_array = []
        x_array = np.array([])
        
        # Select input mode
        if use_query:
            iterator = self.query
        else:
            iterator = self.objects
        
        # Loop
        ##f = file('/nfs/pau/PAUdm/runs/isevilla/analysis/out_data/histo.dat', 'w')
        counter = 1
        counter_found = 1
        if format == None: 
            for db_object in iterator:
                value = db_object.__getattribute__(parameter)
                if parameter == 'mag_auto':
                    image = model.session.query(model.Image).join(model.Detection).filter(model.Detection.id == db_object.id).one()
                    mosaic = model.session.query(model.Mosaic).join(model.Image).filter(model.Image.id == image.id).one()
                    zp_nightly = image.zp_nightly
                    if type(zp_nightly).__name__=='NoneType':
                        continue
                    airmass = mosaic.air_mass
                    value = value + zp_nightly*airmass
                x_array = np.append(x_array,value)
        else:
            x_array = parameter
            #x_array.append(value)
            ##print >> f,value
        ##f.close()
        
        # Create Histogram plot 
        ##x_array = np.array([])
        ##x_array = np.genfromtxt('/nfs/paus/PAUdm/runs/isevilla/analysis/out_data/histo.dat', dtype=None)
        plt.hist(x_array,bins,range=range)        
        
        # Configure plot
        if xlabel == None:
            xlabel = parameter
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.title(title)
        
        # Show or save
        if filename == None:
            plt.show()
        else:
            plt.savefig(filename)
    

            
            
class Phot_Calibrator(object):
    
    def __init__(self, db_mosaic, use_zp_nightly=True):
        """ Load all the zeropoint information from the DB, keep it here for easy access. """
        self.params = {}
        
        # ZP Phot
        self.params['zp_phot'] = [ 0.0 ] # for fictitious CCD #0
        self.params['zp_phot_err'] = [ 0.0 ]
        
        db_zp_phots = model.session.query(model.ZP_Phot).filter(model.ZP_Phot.id == db_mosaic.zp_phot_id).all()
        
        # If we haven't already assigned the correct ZP_phots to this mosaic, then figure it out.
        if( len(db_zp_phots) == 0 ):
            # First, any ZPs for the right filter tray and the right date
            db_zp_phot_query = model.session.query(model.ZP_Phot).filter(model.ZP_Phot.filtertray==db_mosaic.filtertray,model.ZP_Phot.start_date<db_mosaic.date_obs,((model.ZP_Phot.end_date == None) | (model.ZP_Phot.end_date<db_mosaic.date_obs)) )
            # Now check the production info
            db_production = model.session.query(model.Production).filter(model.Production.id==db_mosaic.production_id).one()
            db_zp_phot_query = db_zp_phot_query.join(model.Production)
            db_zp_phot_query = db_zp_phot_query.filter(model.Production.project==db_production.project,model.Production.origin==db_production.origin,model.Production.red_release==db_production.red_release)
            db_zp_phots = db_zp_phot_query.all()
            if( len(db_zp_phots) == 0 ):
                error_msg = "No appropriate ZP_Phot entries found!  Do you need to run the commissioning pipeline?"
                log.error(error_msg)
                raise Exception, error_msg
            if( len(db_zp_phots) > 1 ):
                error_msg = "More than one appropriate ZP Phot entry found!!"
                log.error(error_msg)
                raise Exception, error_msg
            # Set the ZP_Phot ID in the Mosaic
            db_mosaic.zp_phot_id = db_zp_phots[0].id
            
        if len(db_zp_phots) > 1:
            error_msg = "More than one ZP Phot entry has our zp_phot_id of %d !!" %db_mosaic.zp_phot_id
            log.error(error_msg)
            raise Exception, error_msg
            
        for i in range(1,19):
            zpname = "zp_" + "%02d" %i
            zp_err_name = "zp_err_" + "%02d" %i
            self.params['zp_phot'].append( db_zp_phots[0].__getattribute__(zpname) )
            self.params['zp_phot_err'].append( db_zp_phots[0].__getattribute__(zp_err_name) )
            
            
        # ZP Nightly
        self.params['zp_nightly'] = [ 0.0 ] # for fictitious CCD #0
        self.params['zp_nightly_err'] = [ 0.0 ]
        for i in range(1,19):
            if use_zp_nightly:
                db_image = model.session.query(model.Image).filter(model.Image.mosaic==db_mosaic).filter(model.Image.ccd_num==i).one()
                if( db_image.zp_nightly is None ):
                    self.params['zp_nightly'].append( 0.0 )
                    self.params['zp_nightly_err'].append( 0.0 )
                else:
                    self.params['zp_nightly'].append( db_image.zp_nightly )
                    self.params['zp_nightly_err'].append( db_image.zp_nightly_err )
            else:
                self.params['zp_nightly'].append( 0.0 )
                self.params['zp_nightly_err'].append( 0.0 )
                
        # Airmass
        self.params['airmass'] = db_mosaic.airmass
        
        # Softening parameter b if we do asinh mags.....
    
    def flux_to_mag( self, flux, flux_err, ccd_num ):
        (calib_flux, calib_flux_err) = self.calibrate_flux( flux, flux_err, ccd_num )
        mag = 99.0
        mag_error = 99.0
        if( calib_flux > 0.0 ):
            mag = constants.mag_zp-2.5*math.log10(calib_flux)
            mag_error = 1.0857*calib_flux_err/calib_flux
        return (mag, mag_error)


    def calibrate_flux( self, flux, flux_err, ccd_num ):
        zp_phot = self.params['zp_phot'][ccd_num]
        zp_phot_err = self.params['zp_phot_err'][ccd_num]
        zp_nightly = self.params['zp_nightly'][ccd_num]
        zp_nightly_err = self.params['zp_nightly_err'][ccd_num]
        airmass = self.params['airmass']
        
        calib_flux = flux * 10.0**(0.4*(zp_phot + zp_nightly*airmass))
        calib_flux_err = 0.0
        if( flux != 0.0 ):
            calib_flux_err = math.sqrt( (flux_err*flux_err)/(flux*flux) + 0.8483*zp_phot_err*zp_phot_err + 0.8483*zp_nightly_err*zp_nightly_err )*calib_flux
        
        return (calib_flux, calib_flux_err)


class Spectrum(object):
    '''
    Defines a Spectrum object.
    
    Attributes:
      - wavelengths
      - fluxes
      - normalization (i.e., the area under the curve;  useful when this is implemented as a filter curve.)
      
      usage example: 
        spectrum.wavelengths -> numpy array of wavelengths (angstroms)
        spectrum.fluxes -> numpy array of fluxes, same length as wavelengths.
    
    Built-in Methods:
      - flux = filter_curve_spectrum.sed_to_flux( sed_spectrum, norm_filter, norm_flux, single_dw )
        returns the flux of an object with the given SED in this filter.
    '''
    
    def __init__(self, ws=None, fs=None, name=None, file=None):
        self.wavelengths = ws
        self.fluxes = fs
        self.filter_name = name
        self.norm = None
        
        if file is not None:
            #file = open(filename)
            ws = []
            fs = []
            for line in file:
                splitted_line = line.split()
                if( re.match("#", splitted_line[0]) ):
                    continue
                w = float(splitted_line[0])
                f = float(splitted_line[1])
                ws.append(w)
                fs.append(f)
            np_ws = np.array(ws)
            np_fs = np.array(fs)
            
            if self.wavelengths is not None:
                self.fluxes = np.interp( self.wavelengths, np_ws, np_fs )
            else:
                self.wavelengths = np_ws
                self.fluxes = np_fs
                
    
    def rebin( self, wvls ):
        new_fluxes = np.interp( wvls, self.wavelengths, self.fluxes )
        self.wavelengths = wvls.copy()
        self.fluxes = new_fluxes
    
    
    def sed_to_flux( self, sed, norm_filter, norm_mag, single_dw=False ):
        
        norm_result = norm_filter.combine_sed_filter( sed, True )
        
        norm_flux = 10**( (constants.mag_zp-norm_mag)/2.5 )
        
        # find the normalization factor to get the correct flux
        normalization = norm_flux/norm_result
                        
        # multiply the desired filter with the SED and add up the "flux"
        filter_result = self.combine_sed_filter( sed, True )
                                
        # multiply by the normalization factor
        filter_result *= normalization
                                
        return filter_result
        
        
    
    # was combine_sed_filter_matchedwvls in simUtils
    # should be called by a filter curve.
    # sed is a Spectrum type too.
    def combine_sed_filter( self, sed, single_dw=False ):
        
        if( len(sed.wavelengths) != len(self.wavelengths) ):
            error_msg = "Filter and SED arrays are not the same length! len(sed.wavelengths) = %d, len(self.wavelengths) = %d" %(len(sed.wavelengths),len(self.wavelengths))
            log.error(error_msg) 
            
        # basically copying the few equations from 
        # http://arxiv.org/abs/astro-ph/0701508
        
        if self.norm is None:
            error_msg = "Filter normalization is not set!"
            log.error(error_msg)

        if( single_dw ):
            dlambda = self.wavelengths[1] - self.wavelengths[0]
            mult_array = sed.fluxes*sed.wavelengths*self.fluxes # HERE "should" be *wvl, not /wvl, but the quasar is then not blue....
            flux = np.sum(mult_array)*dlambda/self.norm
            
        else:
            dlambda_np_array = np.zeros(len(self.wavelengths))
            dlambda_np_array[0] = self.wavelengths[1] - self.wavelengths[0]
            for i in range(1,len(self.wavelengths)-1):
                dlambda_np_array[i] = (self.wavelengths[i+1]-self.wavelengths[i-1])/2.0 
            dlambda_np_array[len(self.wavelengths)-1] = self.wavelengths[len(self.wavelengths)-1] - self.wavelengths[len(self.wavelengths)-2]
            
            mult_array = sed.fluxes*sed.wavelengths*self.fluxes * dlambda_np_array
            flux = np.sum(mult_array) / self.norm
            
        flux *= 3e18/3631e23
        
        return flux
    


class FilterSystem(object):
    '''
    Defines a standard FilterSystem object.
    
    parameters: system name (e.g. "sdss")
    
    Attributes:
      system -> system name (e.g. "sdss")
      filters -> dictionary of Spectrum objects, with dictionary 
                 indices of the filter names (e.g. 'u')
      
      usage example: 
        filtersystem.filters['g'].wavelength -> list of wavelengths (angstroms)
        
    Built-in Methods:
      - filters = read_filter_curves( system, wvl_list )
        constructs the "filters" dictionary for the given system, 
        with a wavelength array given by wvl_list
    '''
    
    def __init__(self, sysname, instrument, add_efficiency=True):
        self.system = sysname
        
        # a good range, given our input files
        # matches the one in StarSEDs
        #wvl_list = np.arange(95.0,12000.0,2.0) 
        wvl_list = np.arange(95.0,14500.0,2.0) 
        #wvl_list = np.arange(95.0,18500.0,2.0) 
        
        self.filters = self.read_filter_curves( self.system, wvl_list, instrument, add_efficiency )
    
    
    
    def read_filter_curves( self, system, wvl_list, instrument, add_efficiency=True ):
        
        filter_curves = {}
        filter_curves_matchedwvls = {}
        
        atmosphere_file = open(resource_filename('paudm.resources','atmosphere/transmission.sdss.stdunits.dat'))
        qe_file = open(resource_filename('paudm.resources','detector/HPK.real.stdunits.dat'))
        wht_throughput_file = open(resource_filename('paudm.resources','telescope/WHT_throughput_stdunits.dat'))
        
        if( system == "pau" ):
            
            # load up some config info (about the filters)
            for filter_name in instrument['BROAD_FILTERS'] + instrument['NARROW_FILTERS']:
                filter_curves[filter_name] = Spectrum([], [])
                filter_curves_matchedwvls[filter_name] = Spectrum(np.zeros(len(wvl_list)),np.zeros(len(wvl_list)))
            
            for filter_name in instrument['BROAD_FILTERS'] + instrument['NARROW_FILTERS']:
                filter_file = open(resource_filename('paudm.resources', 'filters/pau/'+instrument['FILTER_ASSOCIATIONS'][filter_name]['FILENAME']))
                for line in filter_file:
                    if re.search("(\S+)", line) is None:
                        continue
                    splitted_line = line.split()
                    if re.match("#", splitted_line[0]):
                        continue
                    else:
                        filter_curves[filter_name].wavelengths.append(float(splitted_line[0]))
                        filter_curves[filter_name].fluxes.append(float(splitted_line[1]))
                filter_file.close()
                filter_curves[filter_name].wavelengths = np.array( filter_curves[filter_name].wavelengths )
                filter_curves[filter_name].fluxes = np.array( filter_curves[filter_name].fluxes )
                
                # pad them with some zeroes so that the interpolation doesn't add flux at the edges
                for iter in range(2):
                    last_element = filter_curves[filter_name].wavelengths[len(filter_curves[filter_name].wavelengths)-1]
                    new_list = np.roll( filter_curves[filter_name].wavelengths, 1 )
                    new_list = np.append( new_list, last_element )
                    new_list = np.append( new_list, (last_element-new_list[len(new_list)-2]) + last_element )
                    new_list[0] = new_list[1] - (new_list[2]-new_list[1])
                    filter_curves[filter_name].wavelengths = new_list
                for iter in range(2):
                    last_element = filter_curves[filter_name].fluxes[len(filter_curves[filter_name].fluxes)-1]
                    new_list = np.roll( filter_curves[filter_name].fluxes, 1 )
                    new_list = np.append( new_list, last_element )
                    new_list = np.append( new_list, 0.0 )
                    new_list[0] = 0.0
                    filter_curves[filter_name].fluxes = new_list
                    
                # now make an array with the same wvls as the seds
                filter_curves_matchedwvls[filter_name].wavelengths = wvl_list.copy()
                filter_curves_matchedwvls[filter_name].fluxes = np.interp(wvl_list, filter_curves[filter_name].wavelengths, filter_curves[filter_name].fluxes)
                
            if( add_efficiency ):
                
                # read in the curves
                atmospheric_transmission = Spectrum( ws=wvl_list.copy(), file=atmosphere_file )
                qe_curve = Spectrum( ws=wvl_list.copy(), file=qe_file )
                throughput_curve = Spectrum( ws=wvl_list.copy(), file=wht_throughput_file )
                efficiency = atmospheric_transmission.fluxes*qe_curve.fluxes*throughput_curve.fluxes
                #for i in range( len(throughput_curve.wavelengths) ):
                #    print "%f %f %f %f %f" %(throughput_curve.wavelengths[i], atmospheric_transmission.fluxes[i], qe_curve.fluxes[i], throughput_curve.fluxes[i], efficiency[i])
                
                for filter_name in instrument['BROAD_FILTERS'] + instrument['NARROW_FILTERS']:
                    filter_curves_matchedwvls[filter_name].fluxes *= efficiency
                
        elif( system == "sdss" ):
            sdss_filter_names = ['u', 'g', 'r', 'i', 'z']
            for filter_name in sdss_filter_names:
                filter_curves[filter_name] = Spectrum([],[])
                filter_curves_matchedwvls[filter_name] = Spectrum(np.zeros(len(wvl_list)),np.zeros(len(wvl_list)))
            for filter_name in sdss_filter_names:
                filter_file = open(resource_filename('paudm.resources','filters/sdss/'+ filter_name+'.dat'))
                for line in filter_file:
                    if re.search("(\S+)", line) is None:
                        continue
                    splitted_line = line.split()
                    if re.match("#", splitted_line[0]):
                        continue
                    else:
                        filter_curves[filter_name].wavelengths.append(float(splitted_line[0]))
                        filter_curves[filter_name].fluxes.append(float(splitted_line[1]))
                filter_file.close()
                
                # pad them with some zeroes so that the interpolation doesn't add flux at the edges
                for iter in range(2):
                    last_element = filter_curves[filter_name].wavelengths[len(filter_curves[filter_name].wavelengths)-1]
                    new_list = np.roll( filter_curves[filter_name].wavelengths, 1 )
                    new_list = np.append( new_list, last_element )
                    new_list = np.append( new_list, (last_element-new_list[len(new_list)-2]) + last_element )
                    new_list[0] = new_list[1] - (new_list[2]-new_list[1])
                    filter_curves[filter_name].wavelengths = new_list
                for iter in range(2):
                    last_element = filter_curves[filter_name].fluxes[len(filter_curves[filter_name].fluxes)-1]
                    new_list = np.roll( filter_curves[filter_name].fluxes, 1 )
                    new_list = np.append( new_list, last_element )
                    new_list = np.append( new_list, 0.0 )
                    new_list[0] = 0.0
                    filter_curves[filter_name].fluxes = new_list
                
                # now make an array with the same wvls as the seds
                filter_curves_matchedwvls[filter_name].wavelengths = wvl_list.copy()
                filter_curves_matchedwvls[filter_name].fluxes = np.interp(wvl_list, filter_curves[filter_name].wavelengths, filter_curves[filter_name].fluxes)
                
        
        elif( system == "johnson" ):
            johnson_filter_names = ['I']
            for filter_name in johnson_filter_names:
                filter_curves[filter_name] = Spectrum([],[])
                filter_curves_matchedwvls[filter_name] = Spectrum(np.zeros(len(wvl_list)),np.zeros(len(wvl_list)))
            for filter_name in johnson_filter_names:
                filter_file = open(resource_filename("paudm.resources","filters/johnson/"+filter_name + ".pb"))
                for line in filter_file:
                    if re.search("(\S+)", line) is None:
                        continue
                    splitted_line = line.split()
                    if re.match("#", splitted_line[0]):
                        continue
                    else:
                        filter_curves[filter_name].wavelengths.append(float(splitted_line[0]))
                        filter_curves[filter_name].fluxes.append(float(splitted_line[1]))
                filter_file.close()
                
                # pad them with some zeroes so that the interpolation doesn't add flux at the edges
                for iter in range(2):
                    last_element = filter_curves[filter_name].wavelengths[len(filter_curves[filter_name].wavelengths)-1]
                    new_list = np.roll( filter_curves[filter_name].wavelengths, 1 )
                    new_list = np.append( new_list, last_element )
                    new_list = np.append( new_list, (last_element-new_list[len(new_list)-2]) + last_element )
                    new_list[0] = new_list[1] - (new_list[2]-new_list[1])
                    filter_curves[filter_name].wavelengths = new_list
                for iter in range(2):
                    last_element = filter_curves[filter_name].fluxes[len(filter_curves[filter_name].fluxes)-1]
                    new_list = np.roll( filter_curves[filter_name].fluxes, 1 )
                    new_list = np.append( new_list, last_element )
                    new_list = np.append( new_list, 0.0 )
                    new_list[0] = 0.0
                    filter_curves[filter_name].fluxes = new_list
                    
                # now make an array with the same wvls as the seds
                filter_curves_matchedwvls[filter_name].wavelengths = wvl_list.copy()
                filter_curves_matchedwvls[filter_name].fluxes = np.interp(wvl_list, filter_curves[filter_name].wavelengths, filter_curves[filter_name].fluxes)
                
            if( add_efficiency ):
                
                # read in the curves
                atmospheric_transmission = Spectrum( ws=wvl_list.copy(), file=atmosphere_file )
                qe_curve = Spectrum( ws=wvl_list.copy(), file=qe_file )
                throughput_curve = Spectrum( ws=wvl_list.copy(), file=wht_throughput_file )
                efficiency = atmospheric_transmission.fluxes*qe_curve.fluxes*throughput_curve.fluxes
                #for i in range( len(throughput_curve.wavelengths) ):
                #    print "%f %f %f %f %f" %(throughput_curve.wavelengths[i], atmospheric_transmission.fluxes[i], qe_curve.fluxes[i], throughput_curve.fluxes[i], efficiency[i])
                
                for filter_name in johnson_filter_names:
                    filter_curves_matchedwvls[filter_name].fluxes *= efficiency
                    
                
        elif( system == "bessell" ):
            #bessell_filter_names = ['Ux', 'B', 'V', 'Rc', 'I', 'J', 'H', 'K']
            #bessell_filter_names = ['Ux', 'B', 'V', 'Rc', 'Ic', 'J', 'H']
            bessell_filter_names = ['Ux', 'B', 'V', 'Rc', 'Ic', 'J']
            for filter_name in bessell_filter_names:
                filter_curves[filter_name] = Spectrum([],[])
                filter_curves_matchedwvls[filter_name] = Spectrum(np.zeros(len(wvl_list)),np.zeros(len(wvl_list)))
            for filter_name in bessell_filter_names:
                filter_file = open(resource_filename("paudm.resources", "filters/bessell/" +filter_name + ".pb"))
                for line in filter_file:
                    if re.search("(\S+)", line) is None:
                        continue
                    splitted_line = line.split()
                    if re.match("#", splitted_line[0]):
                        continue
                    else:
                        filter_curves[filter_name].wavelengths.append(float(splitted_line[0]))
                        filter_curves[filter_name].fluxes.append(float(splitted_line[1]))
                filter_file.close()
                
                # pad them with some zeroes so that the interpolation doesn't add flux at the edges
                for iter in range(2):
                    last_element = filter_curves[filter_name].wavelengths[len(filter_curves[filter_name].wavelengths)-1]
                    new_list = np.roll( filter_curves[filter_name].wavelengths, 1 )
                    new_list = np.append( new_list, last_element )
                    new_list = np.append( new_list, (last_element-new_list[len(new_list)-2]) + last_element )
                    new_list[0] = new_list[1] - (new_list[2]-new_list[1])
                    filter_curves[filter_name].wavelengths = new_list
                for iter in range(2):
                    last_element = filter_curves[filter_name].fluxes[len(filter_curves[filter_name].fluxes)-1]
                    new_list = np.roll( filter_curves[filter_name].fluxes, 1 )
                    new_list = np.append( new_list, last_element )
                    new_list = np.append( new_list, 0.0 )
                    new_list[0] = 0.0
                    filter_curves[filter_name].fluxes = new_list
                    
                # now make an array with the same wvls as the seds
                filter_curves_matchedwvls[filter_name].wavelengths = wvl_list.copy()
                filter_curves_matchedwvls[filter_name].fluxes = np.interp(wvl_list, filter_curves[filter_name].wavelengths, filter_curves[filter_name].fluxes)
                
            if( add_efficiency ):
                
                # read in the curves
                atmospheric_transmission = Spectrum( ws=wvl_list.copy(), file=atmosphere_file )
                qe_curve = Spectrum( ws=wvl_list.copy(), file=qe_file )
                throughput_curve = Spectrum( ws=wvl_list.copy(), file=wht_throughput_file )
                efficiency = atmospheric_transmission.fluxes*qe_curve.fluxes*throughput_curve.fluxes
                #for i in range( len(throughput_curve.wavelengths) ):
                #    print "%f %f %f %f %f" %(throughput_curve.wavelengths[i], atmospheric_transmission.fluxes[i], qe_curve.fluxes[i], throughput_curve.fluxes[i], efficiency[i])
                
                for filter_name in bessell_filter_names:
                    filter_curves_matchedwvls[filter_name].fluxes *= efficiency
                    
        # this includes the atmosphere, detector QE, and WHT throughput.
        elif( system == "pau_efficiency" ):
            
            # read in the curves
            atmospheric_transmission = Spectrum( ws=wvl_list.copy(), file=atmosphere_file )
            qe_curve = Spectrum( ws=wvl_list.copy(), file=qe_file )
            throughput_curve = Spectrum( ws=wvl_list.copy(), file=wht_throughput_file )
            
            filter_curves_matchedwvls['qe'] = Spectrum( ws=wvl_list.copy(), fs=qe_curve.fluxes.copy() )
            filter_curves_matchedwvls['throughput'] = Spectrum( ws=wvl_list.copy(), fs=throughput_curve.fluxes.copy() )
            filter_curves_matchedwvls['total'] = Spectrum( ws=wvl_list.copy(), fs=atmospheric_transmission.fluxes*qe_curve.fluxes*throughput_curve.fluxes )
            
            
        else:
            error_msg = "pau, sdss, johnson, and bessell filter systems, and pau_efficiency, are implemented for now, not %s!" %system
            log.error(error_msg)
            
            
        # calculate the normalization for each filter, 
        # which we will need for making mags out of seds+filters
        for filter_name in filter_curves_matchedwvls.keys():
            dlambda = filter_curves_matchedwvls[filter_name].wavelengths[1] - filter_curves_matchedwvls[filter_name].wavelengths[0]
            mult_array = filter_curves_matchedwvls[filter_name].fluxes/filter_curves_matchedwvls[filter_name].wavelengths
            filter_curves_matchedwvls[filter_name].norm = np.sum(mult_array)*dlambda
        
        atmosphere_file.close()
        qe_file.close()
        wht_throughput_file.close()    
        return filter_curves_matchedwvls
    


class StarSEDs(object):
    '''
    Defines an object that holds template stellar SEDs.
    
    parameters: none
    
    Attributes:
      - seds -> a dictionary of Spectrum objects that are the stellar SEDs
      
      usage example: 
      
      
    Built-in Methods:
      - interpolate_seds -> interpolates between the Pickles templates
    '''
    
    def __init__(self):
        self.seds = dict()
        
        sed_files = resource_listdir('paudm.resources','seds/stars')
        
        # matches the one in FilterSystem
        #wvl_list = np.arange(95.0,12000.0,2.0)
        wvl_list = np.arange(95.0,14500.0,2.0)
        #wvl_list = np.arange(95.0,18500.0,2.0)
        re_sed = re.compile("(\S+)")
        re_comm = re.compile("#")
        for index, file in enumerate(sed_files):
            if file[-3:] != "sed":
                continue
            self.seds[file] = Spectrum([],[])
            
            sed_file = open(resource_filename ('paudm.resources','seds/stars/'+file))
            for line in sed_file:
                #if re.search("(\S+)", line) is None:
                if not re_sed.search(line):
                    continue
                splitted_line = line.split()
                if re_comm.match(splitted_line[0]):
                    continue
                else:
                    self.seds[file].wavelengths.append(float(splitted_line[0]))
                    self.seds[file].fluxes.append(float(splitted_line[1]))
            sed_file.close()
            # a few SEDs are not what we want, e.g. zodiacal, flat
            # these are a different length from the others.  remove them.
            if len(self.seds[file].wavelengths) != 4773:
                del self.seds[file]
            else:
                self.seds[file].fluxes = np.interp( wvl_list, self.seds[file].wavelengths, self.seds[file].fluxes )
                self.seds[file].wavelengths = wvl_list 
            
        log.debug("PAUStars: read in Pickles SEDs.")
        log.debug("Interpolating SEDs... starting from %d SEDs." %len(self.seds))
        self.interpolateSEDs()
        log.debug("Done interpolating.  Now %d SEDs." %len(self.seds))
    
    
    def interpolateSEDs( self ):
        luminosities = ( 'i', 'ii', 'iii', 'iv', 'v' )
        letters = ( 'o', 'b', 'a', 'f', 'g', 'k', 'm' )
        numbers = ( '0', '1', '2', '3', '4', '5', '6', '7', '8', '9' )
        
        oldsize = 0
        newsize = 1
        iterations = 0
        while oldsize != newsize:
            
            oldsize = newsize
            
            for lu in range(len(luminosities)):
                luminosity = luminosities[lu]
            # first interpolate in the letters+numbers direction
                start_le = None
                start_nu = None
                le = 0
                while le < len(letters):
                    letter = letters[le]
                    nu = 0
                    while nu < len(numbers):
                        number = numbers[nu]
                        sedName = letter + number + luminosity + ".sed"
                        if sedName in self.seds:
                            if start_le is None:
                                start_le = le
                                start_nu = nu
                            else:
                                if start_le == le and start_nu == nu-1:
                                    # we don't need to interpolate
                                    start_nu = nu
                                elif start_le == le-1 and start_nu == 9 and nu == 0:
                                    # we don't need to interpolate
                                    start_le = le
                                    start_nu = nu
                                else:
                                    # we do need to interpolate!
                                    # the previous one has the same letter
                                    if start_le == le:
                                        for step_back in range(1,nu-start_nu):
                                            step_forward = nu-start_nu - step_back
                                            weight_back = (1./float(step_back)) / (1./float(step_back) + 1./float(step_forward))
                                            weight_forward = (1./float(step_forward)) / (1./float(step_back)+1./float(step_forward))
                                            name_back = letters[start_le] + numbers[start_nu] + luminosity + ".sed"
                                            new_name = letters[start_le] + numbers[start_nu+step_back] + luminosity + ".sed"
                                            self.seds[new_name] = self.combine_seds( self.seds[name_back],weight_back,self.seds[sedName],weight_forward)
                                            
                                        start_le = le
                                        start_nu = nu
                                    # if the previous one has a previous letter
                                    else:
                                        n_steps = 9-start_nu + nu + 1
                                        for step_back in range(1,n_steps):
                                            step_forward = n_steps - step_back
                                            weight_back = (1./float(step_back)) / (1./float(step_back) + 1./float(step_forward))
                                            weight_forward = (1./float(step_forward)) / (1./float(step_back)+1./float(step_forward))
                                            name_back = letters[start_le] + numbers[start_nu] + luminosity + ".sed"
                                            if start_nu+step_back > 9:
                                                new_name = letters[le] + numbers[nu-step_forward] + luminosity + ".sed"
                                            else:
                                                new_name = letters[start_le] + numbers[start_nu+step_back] + luminosity + ".sed"
                                            self.seds[new_name] = self.combine_seds(self.seds[name_back],weight_back,self.seds[sedName],weight_forward)
                                        start_le = le
                                        start_nu = nu
                        else:
                            pass
                        
                        nu = nu+1
                    le = le+1
                    
        # now interpolate in the i ii iii iv v direction.
            for le in range(len(letters)):
                letter = letters[le]
                for nu in range(len(numbers)):
                    number = numbers[nu]
                    start_lu = None
                    for lu in range(len(luminosities)):
                        luminosity = luminosities[lu]
                        sedName = letter + number + luminosity + ".sed"
                        if sedName in self.seds:
                            if start_lu is None:
                                start_lu = lu
                            else:
                                if start_lu == lu-1:
                                    # we don't need to interpolate
                                    start_lu = lu
                                else:
                                    # we do need to interpolate!
                                    for step_back in range(1,lu-start_lu):
                                        step_forward = lu-start_lu - step_back
                                        weight_back = (1./float(step_back)) / (1./float(step_back) + 1./float(step_forward))
                                        weight_forward = (1./float(step_forward)) / (1./float(step_back)+1./float(step_forward))
                                        name_back = letter + number + luminosities[start_lu] + ".sed"
                                        new_name = letter + number + luminosities[start_lu+step_back] + ".sed"
                                        self.seds[new_name] = self.combine_seds( self.seds[name_back],weight_back,self.seds[sedName],weight_forward)
                                        
                                    start_lu = lu
            newsize = len(self.seds)
            iterations = iterations+1
            
        log.debug("Finished after %d iterations" %iterations)
        
    
    def combine_seds( self, sed1, weight1, sed2, weight2 ):
        if len(sed1.wavelengths) != len(sed2.wavelengths):
            error_msg = "combine_seds: SEDs are different lengths!"
            log.error(error_msg)
            raise Exception, error_msg
            
        newSED = Spectrum(np.zeros(len(sed1.wavelengths)),np.zeros(len(sed1.wavelengths)))
        newSED.wavelengths = sed1.wavelengths
        newSED.fluxes = weight1*sed1.fluxes + weight2*sed2.fluxes
        
        return newSED
    


