#!/usr/bin/env python
# -*- coding: utf-8 -*-

import hashlib
import optparse

from sqlalchemy import ForeignKeyConstraint, Index, PrimaryKeyConstraint, UniqueConstraint
from sqlalchemy import create_engine
from sqlalchemy.types import BigInteger, Boolean, Date, Enum, Float, Integer, SmallInteger, String, Text, Time, DateTime
from sqlalchemy.ext.compiler import compiles
from sqlalchemy.orm import mapper, relationship, sessionmaker, scoped_session, synonym, contains_eager
from zope.sqlalchemy import ZopeTransactionExtension
from sqlalchemy.orm.collections import attribute_mapped_collection
from sqlalchemy.orm.exc import NoResultFound, MultipleResultsFound #@UnusedImport
from sqlalchemy.orm import backref
from sqlalchemy.sql.expression import func
from sqlalchemy.ext.declarative import declarative_base
from base import Column, MetaData, Table
from sqlalchemy.ext.hybrid import hybrid_property
import brownthrower.model 

Base = brownthrower.model.Base
metadata      = Base.metadata
session       = None
tables        = {}
mappers       = {}

_salt = 'djhferuyunniurnlp097jlknf8holanhgnhhrf'

def init(url):
    """
    Initialize the tables, mapping classes and establish a session with the DB
    
    @param url: Configuration URL for the SQLAlchemy engine
    @type url: str
    """
    global session
    

    twophase = True    
    if url.startswith('postgresql'):
        engine = create_engine(url, echo = False)
    else:
        engine = create_engine(url, connect_args={'check_same_thread':False}, echo = False)
        twophase = False
    metadata.bind = engine
    
    session = scoped_session(sessionmaker(bind=engine,
        twophase = twophase,
       #extension = ZopeTransactionExtension()
       ))()
    
    return session
'''
class User(Base):
    __tablename__ = 'user'
    __table_args__ = (
        # Primary key
        PrimaryKeyConstraint('id'),
        # Unique key
        UniqueConstraint('email'),
    )
    # Columns
    id        = Column(Integer,     nullable=False )
    email     = Column(String(32),  nullable=False,  index=True)
    name      = Column(String(64),  nullable=True)
    surname   = Column(String(64),  nullable=False  )
    _password = Column(String(128), nullable=False, name='password')
    permissions   = Column(Integer,     nullable=False)
    validated  = Column(Boolean,     nullable=False)
    password = synonym('_password', map_column=True)
    #'password' : synonym('_password', map_column=True),     
#    @synonym_for("_password")
#    def name(self):
#        return "password: %s" % _password    
    def __init__(self, email, name, surname, password, permissions):
        self.email = email
        self.name = name
        self.surname = surname 
        self._password = hashlib.sha512(password + _salt).hexdigest()
        self.validated = False
        self.permissions = permissions
        
    def _get_password(self):
        return self._password

    def _set_password(self, password):
        self._password = hashlib.sha512(password + _salt).hexdigest()
    
    password = property(_get_password, _set_password)
    
    @classmethod
    def get_by_username(cls, username):
        return session.query(cls).filter(cls.email == username).first()
    
    @classmethod
    def check_password(cls, username, password):
        user = cls.get_by_username(username)
        if not user:
            return False 
        hash_password = hashlib.sha512(password + _salt).hexdigest()
        return  hash_password == user.password
    
    @classmethod
    def check_validation(cls, username):
        user = cls.get_by_username(username)
        if not user:
            return False
        if user.validated == 1:
            return True 
        return False
    @hybrid_property
    def password(self):
        return self._password
'''
class Obs_set(Base):
    __tablename__ = 'obs_set'
    __table_args__ = (
        # Primary key
        PrimaryKeyConstraint('id'),
        # Unique key
        UniqueConstraint('obs_set', 'instrument'),
    )
    # Columns
    id        = Column(Integer,     nullable=False )   # 
    instrument = Column(String(16),  nullable=False)    # Name of the instrument used at the observation (i.e. PAUCam1.0)
    log =  Column(String(128),  nullable=True)        # Camera XML log file path
    rjd_start = Column(Float(53),  nullable=False, index = True)     # First exposure observation time from current obs_set
    rjd_stop = Column(Float(53),  nullable=False)     #Last exposure observation time from current obs_set
    night = Column(Date,  nullable=False)             #Observation date
    notes = Column(Text,  nullable=False)             #Observer's notes and comments
    operator = Column(String(128),  nullable=False)    #Observer's name
    obs_set = Column(String(128),  nullable=False)      # Observation Set identifier, from header
    
    #Relationships

    mosaics      = relationship('Mosaic',           back_populates="obs_set")
   
    #The obs_set table tracks observation nights. Each entry corresponds to a night of observation

class Production(Base): 
    __tablename__ = 'production'
    __table_args__ = (
        # Primary key
        PrimaryKeyConstraint('id'),
        # Unique key
        UniqueConstraint('pipeline', 'release'),
        ForeignKeyConstraint(['input_production_id'], ['production.id'], onupdate='CASCADE', ondelete='CASCADE'),
    )
    # Columns
    id                  = Column(Integer,    nullable=False ) #
    input_production_id = Column(Integer,    nullable=True )  # input release id (from configuration)
    pipeline            = Column(String(16), nullable=True)   # Pipeline Name [pixelsim, nightly, memba, analysis]
    release             = Column(String(16), nullable=False)  # Major release name [TESTX, DRX]
    software_version    = Column(String(16), nullable=False)  # Package version [DCX, vX.X]
    job_id              = Column(Integer,    nullable=False)  # id of the job that generates the current production

    #Relationships
    mosaics        = relationship('Mosaic',         back_populates="production")
    zp_phots       = relationship('Zp_phot',        back_populates="production")
    global_objects = relationship('Global_object',  back_populates="production")
    coadd_objects  = relationship('Coadd_object',   back_populates="production")
    photo_zs       = relationship('Photo_z',        back_populates="production")
    truth_objects  = relationship('Truth_Object',   back_populates="production")
    targets        = relationship('Target',         back_populates="production")
    parent         = relationship('Production',     back_populates = 'children',
                         primaryjoin = 'Production.input_production_id == Production.id',
                         remote_side = 'Production.id')
    children       = relationship('Production',     back_populates = 'parent',
                         primaryjoin = 'Production.input_production_id == Production.id',
                         cascade     = 'all, delete-orphan', passive_deletes = True)
    job            = relationship('Job',
                          primaryjoin = 'Production.job_id == Job.id',
                          foreign_keys = '[ Production.job_id ]',
                          remote_side = '[ Job.id ]',
                          uselist = False,
                          backref = backref("production", uselist = False))
    
    
    #The production table holds information about the hardware and software used to process the data.\n"
    #"Each time the software is updated, configuration changes or computing system has been modified, a new production entry is stored.   
    
    
class Zp_phot(Base): 
    __tablename__ = 'zp_phot'
    __table_args__ = (
        # Primary key
        PrimaryKeyConstraint('id'),
        # Unique key
        UniqueConstraint('production_id', 'id'),
        ForeignKeyConstraint(['production_id'], ['production.id'], onupdate='CASCADE', ondelete='CASCADE'),
    )
    # Columns
    id        = Column(BigInteger,     nullable=False )   # 
    production_id = Column(Integer,  nullable=False)    # Production number
    start_date =  Column(Date,  nullable=False)        # Valid start date
    end_date = Column(Date,  nullable=True)     # Valid end date
    filtertray = Column(String(16),  nullable=True)     #Filter Tray name
    zp_01 = Column(Float(24),  nullable=False)             #Photometric Zeropoint magnitude (ZPphot) for CCD 01
    zp_02 = Column(Float(24),  nullable=False)             #Photometric Zeropoint magnitude (ZPphot) for CCD 02
    zp_03 = Column(Float(24),  nullable=False)             #Photometric Zeropoint magnitude (ZPphot) for CCD 03
    zp_04 = Column(Float(24),  nullable=False)             #Photometric Zeropoint magnitude (ZPphot) for CCD 04
    zp_05 = Column(Float(24),  nullable=False)             #Photometric Zeropoint magnitude (ZPphot) for CCD 05
    zp_06 = Column(Float(24),  nullable=False)             #Photometric Zeropoint magnitude (ZPphot) for CCD 06
    zp_07 = Column(Float(24),  nullable=False)             #Photometric Zeropoint magnitude (ZPphot) for CCD 07
    zp_08 = Column(Float(24),  nullable=False)             #Photometric Zeropoint magnitude (ZPphot) for CCD 08
    zp_09 = Column(Float(24),  nullable=False)             #Photometric Zeropoint magnitude (ZPphot) for CCD 09
    zp_10 = Column(Float(24),  nullable=False)             #Photometric Zeropoint magnitude (ZPphot) for CCD 10
    zp_11 = Column(Float(24),  nullable=False)             #Photometric Zeropoint magnitude (ZPphot) for CCD 11
    zp_12 = Column(Float(24),  nullable=False)             #Photometric Zeropoint magnitude (ZPphot) for CCD 12
    zp_13 = Column(Float(24),  nullable=False)             #Photometric Zeropoint magnitude (ZPphot) for CCD 13
    zp_14 = Column(Float(24),  nullable=False)             #Photometric Zeropoint magnitude (ZPphot) for CCD 14
    zp_15 = Column(Float(24),  nullable=False)             #Photometric Zeropoint magnitude (ZPphot) for CCD 15
    zp_16 = Column(Float(24),  nullable=False)             #Photometric Zeropoint magnitude (ZPphot) for CCD 16
    zp_17 = Column(Float(24),  nullable=False)             #Photometric Zeropoint magnitude (ZPphot) for CCD 17
    zp_18 = Column(Float(24),  nullable=False)             #Photometric Zeropoint magnitude (ZPphot) for CCD 18
    zp_err_01 = Column(Float(24),  nullable=False)             #Photometric Zeropoint magnitude error (σZPphot)  for CCD 01
    zp_err_02 = Column(Float(24),  nullable=False)             #Photometric Zeropoint magnitude error (σZPphot)  for CCD 02
    zp_err_03 = Column(Float(24),  nullable=False)             #Photometric Zeropoint magnitude error (σZPphot)  for CCD 03
    zp_err_04 = Column(Float(24),  nullable=False)             #Photometric Zeropoint magnitude error (σZPphot)  for CCD 04
    zp_err_05 = Column(Float(24),  nullable=False)             #Photometric Zeropoint magnitude error (σZPphot)  for CCD 05
    zp_err_06 = Column(Float(24),  nullable=False)             #Photometric Zeropoint magnitude error (σZPphot)  for CCD 06
    zp_err_07 = Column(Float(24),  nullable=False)             #Photometric Zeropoint magnitude error (σZPphot)  for CCD 07
    zp_err_08 = Column(Float(24),  nullable=False)             #Photometric Zeropoint magnitude error (σZPphot)  for CCD 08
    zp_err_09 = Column(Float(24),  nullable=False)             #Photometric Zeropoint magnitude error (σZPphot)  for CCD 09
    zp_err_10 = Column(Float(24),  nullable=False)             #Photometric Zeropoint magnitude error (σZPphot)  for CCD 10
    zp_err_11 = Column(Float(24),  nullable=False)             #Photometric Zeropoint magnitude error (σZPphot)  for CCD 11
    zp_err_12 = Column(Float(24),  nullable=False)             #Photometric Zeropoint magnitude error (σZPphot)  for CCD 12
    zp_err_13 = Column(Float(24),  nullable=False)             #Photometric Zeropoint magnitude error (σZPphot)  for CCD 13
    zp_err_14 = Column(Float(24),  nullable=False)             #Photometric Zeropoint magnitude error (σZPphot)  for CCD 14
    zp_err_15 = Column(Float(24),  nullable=False)             #Photometric Zeropoint magnitude error (σZPphot)  for CCD 15
    zp_err_16 = Column(Float(24),  nullable=False)             #Photometric Zeropoint magnitude error (σZPphot)  for CCD 16
    zp_err_17 = Column(Float(24),  nullable=False)             #Photometric Zeropoint magnitude error (σZPphot)  for CCD 17
    zp_err_18 = Column(Float(24),  nullable=False)             #Photometric Zeropoint magnitude error (σZPphot)  for CCD 18  
    
    #Relationships

    production    = relationship('Production',       back_populates="zp_phots")
Index('uk_production_zp_phot', Zp_phot.production_id, Zp_phot.id)

#        "Contains information for the photometric calibration.\n"
#        "Holds a filter tray set (18 CCDs) of ZPphot for an interval of time."
class Mosaic(Base): 
    """
    Defines a standard Mosaic object.
    
    Keys:
     - id
     - production_id
     - obs_set_id
     - zp_phot_id
     
     Fields:
     - filename
     - archivepath
     - kind
     - exp_num
     - obs_title
     - ra
     - dec
     - equinox
     - date_obs
     - time_obs
     - rjd_obs
     - date_creat
     - time_creat
     - exp_time
     - air_mass
     - telfocus
     - instrument
     - filtertray
     - filtertray_tmp
     - nextend
     - guide_enabled
     - guide_fwhm
     - guide_var
     - extinction
     - extinction_err
     - wind_spd
     - wind_dir
     - amb_temp
     - humidity
     - pressure
     - eqa_1
     - eqa_2
     - eqa_3
     - eqa_4
     - eqa_5
    
    Relationships:
     - production
     - obs_set
     - images
    """
    #TODO: evaluate and remove maybe?
#    @ensure_criteria({
#        'mosaic' : ['kind'],
#    })

    __tablename__ = 'mosaic'
    __table_args__ = (
        # Primary key
        PrimaryKeyConstraint('id'),
        # Unique key


        # Foreign key
        ForeignKeyConstraint(['production_id'], ['production.id'],                                     onupdate='CASCADE', ondelete='CASCADE'),
        ForeignKeyConstraint(['obs_set_id'],        ['obs_set.id'],                                            onupdate='CASCADE', ondelete='CASCADE'),
        ForeignKeyConstraint(['production_id', 'zp_phot_id'], ['zp_phot.production_id', 'zp_phot.id'], onupdate='CASCADE', ondelete='RESTRICT'),
        UniqueConstraint('archivepath', 'filename'),
        UniqueConstraint('production_id', 'obs_set_id', 'kind', 'exp_num'),

    )
    # Columns
    id        = Column(BigInteger,     nullable=False )   # 
    production_id = Column(Integer,  nullable=False)    # Production identifier
    obs_set_id =  Column(Integer,  nullable=False)        # obs_set number
    zp_phot_id = Column(Integer,  nullable=True)     # ZPphot identifier
    filename = Column(String(128),  nullable=False)     #File name
    archivepath = Column(String(128),  nullable=False)             #Path in the archive
    kind = Column(Enum('ARC', 'BIAS', 'DARK', 'FLASH', 'FLAT', 'FOCUS', 'GLANCE', 'PUPIL',
                       'SCRATCH', 'SKY', 'TARGET','MBIAS','MFLAT','RED_SCI','RED_WEIGHT','RED_MASK',
                       name='mosaic_kind'),  nullable=False)             #Mosaic image type
    exp_num = Column(Integer,  nullable=False)     # Exposure number
    obs_title = Column(String(32),   nullable=False)     #Observation title
    ra = Column(Float(53),    nullable=False)      #Telescope Right Ascension pointing (deg)
    dec = Column(Float(53),    nullable=False)     # Telescope Declination pointing (deg)
    equinox = Column(String(32),  nullable=False)     #Equinox of telescope coordinates
    date_obs = Column(Date,  nullable=False)  #Observation date
    time_obs = Column(Time,  nullable=False)     # Observation time (Universal time)
    rjd_obs = Column(Float(53),  nullable=False)     #Observation Reduced Modified Julian Day
    date_creat = Column(Date,  nullable=False)      # File creation date
    time_creat = Column(Time,  nullable=False)     #File creation time
    exp_time = Column(Float(24),  nullable=False)     #Exposure time (s)
    airmass = Column(Float(24),  nullable=True)      #airmass
    telfocus = Column(Float(24),  nullable=True)     # Telescope Focus
    instrument = Column(String(32),  nullable=False)     #Instrument name
    filtertray = Column(String(16),  nullable=True)      #Filter Tray name
    filtertray_tmp = Column(Float(24),  nullable=True)     # Filter Tray temperature (deg C)
    nextend = Column(SmallInteger,  nullable=False)     #Filter Tray name
    guide_enabled = Column(Boolean,  nullable=True)   #Number of extensions
    guide_fwhm = Column(Float(24),    nullable=True )     # Seeing FWHM measured at guiding star (arcsec)
    guide_var = Column(Float(24),    nullable=True )     #Flux variance measured at guiding star (counts)
    extinction = Column(Float(24),    nullable=True )      #Extinction value measured in reduction for photometric calibration. Available in reduced image only.
    extinction_err = Column(Float(24),    nullable=True )     # Extinction error value. Available in reduced image only.
    wind_spd = Column(Float(24),    nullable=True )     #Wind speed (kph)
    wind_dir = Column(Float(24),    nullable=True )   #Wind direction (deg)
    amb_temp = Column(Float(24),    nullable=True)     # Ambient temperature (deg C)
    humidity = Column(Float(24),    nullable=True)     #Ambient relative humidity (percent)
    pressure = Column(Float(24),    nullable=True)  #Barometic pressure (mbar)
    eqa_1 = Column(Float(24),    nullable=True)     # Quality Analysis check at observatory 1
    eqa_2 = Column(Float(24),    nullable=True)     # Quality Analysis check at observatory 2
    eqa_3 = Column(Float(24),    nullable=True)      # Quality Analysis check at observatory 3
    eqa_4 = Column(Float(24),    nullable=True)     # Quality Analysis check at observatory 4
    eqa_5 = Column(Float(24),    nullable=True)    # Quality Analysis check at observatory 5
    # Relationships

    production = relationship('Production', back_populates="mosaics")
    obs_set        = relationship('Obs_set',        back_populates="mosaics")
    images = relationship('Image',      back_populates="mosaic")
    @classmethod
    def criteria_query(cls, criteria={}):
        """
        obs_sets a query to retrieve Mosaic instances to be prestaged
        """
        q = session.query(cls).join(cls.production).join(cls.obs_set)
        q = add_criteria(q, criteria)
        return q
Index('ik_location', Mosaic.production_id, Mosaic.ra, Mosaic.dec)
Index('ik_rjdobs', Mosaic.production_id, Mosaic.rjd_obs)
#   "Contains detailed header information for each mosaic.\n"
#        "Each mosaic entry will have several ccd images associated (18 in the case of PAU full focal plane)."


class Image(Base):
    """
    Defines a standard Image object.
    
    Keys:
     - id
     - mosaic_id
    
    Fields:
     - image_num
     - ccd_num
     - amp_num
     - filter
     - wavelength
     - waveband
     - gain
     - rdnoise
     - naxis1
     - naxis2
     - zp_nightly
     - zp_nightly_sigma
     - psf_fwhm
     - cqa_1
     - cqa_2
     - cqa_3
     - cqa_4
     - cqa_5
    
    Relationships:
     - mosaic
     - detections
    """
    __tablename__ = 'image'
    __table_args__ = (
        PrimaryKeyConstraint('id'),
        ForeignKeyConstraint(['mosaic_id'], ['mosaic.id'], onupdate='CASCADE', ondelete='CASCADE'),
        UniqueConstraint('mosaic_id', 'ccd_num', 'amp_num'),
        UniqueConstraint('mosaic_id', 'image_num'),
        )
    # Keys
    id = Column(BigInteger,   nullable=False)
    mosaic_id  = Column( BigInteger,   nullable=False)
    # Fields
    image_num = Column(SmallInteger, nullable=False) #Extension number
    ccd_num = Column(          SmallInteger, nullable=False) #CCD Number"
    amp_num = Column(           SmallInteger, nullable=False) #Amplifier number (-1 for full CCD)
    filter  = Column(      String(4),    nullable=True)  #Filter name"),
    wavelength = Column(Float(24),    nullable=False) #Wavelength at filter center (nm)"),
    waveband = Column(          Float(24),    nullable=False) #
    gain =  Column(              Float(24),    nullable=False) #Detector gain at amplifier (e-/ADU)"),
    rdnoise = Column(           Float(24),    nullable=False) #Amplifier Readout Noise"),
    naxis1 = Column(            SmallInteger, nullable=False) #Axis 1 size (pix)"),
    naxis2= Column(           SmallInteger, nullable=False)  #Axis 2 size (pix)"),
    ra_min = Column(            Float(53),    nullable=True) #Image corner ra min (deg)"),
    ra_max = Column(            Float(53),    nullable=True) #Image corner ra max (deg)"),
    dec_min = Column(           Float(53),    nullable=True) #Image corner dec min (deg)"),
    dec_max = Column(           Float(53),    nullable=True) #Image corner dec min (deg)"),
    zp_nightly = Column(        Float(24),    nullable=True) #Zeropoint magnitude computed at the Nightly pipeline"),
    zp_nightly_err = Column(    Float(24),    nullable=True) #Zeropoint magnitude error computed at the Nightly pipeline"),
    psf_fwhm = Column(          Float(24),    nullable=True) #PSF FWHM measured on image. Available in reduced image only."),
    cqa_1 = Column(             Float(24),    nullable=False) #Quality Analysis check at observatory 1"),
    cqa_2 = Column(             Float(24),    nullable=False) #Quality Analysis check at observatory 2"),
    cqa_3 = Column(             Float(24),    nullable=False) #Quality Analysis check at observatory 3"),
    cqa_4 = Column(             Float(24),    nullable=False) #Quality Analysis check at observatory 4"),
    cqa_5 = Column(             Float(24),    nullable=False) #Quality Analysis check at observatory 5"),
    #Relationships

    mosaic          = relationship('Mosaic',           back_populates="images")
    detections     = relationship('Detection',        back_populates="image")
    detections_by_id = relationship('Detection',       collection_class=attribute_mapped_collection('id'))


    #Contains detailed header information for a single CCD image
Index('ik_imagelocation', Image.ra_min, Image.ra_max, Image.dec_min, Image.dec_max)


class Detection(Base):
    __tablename__ = 'detection'
    __table_args__ = (
        # Constraints
        PrimaryKeyConstraint('id'),
        ForeignKeyConstraint(['image_id'],          ['image.id'],      onupdate='CASCADE', ondelete='RESTRICT'),
        #ForeignKeyConstraint(['global_object_id'],  ['global_object.id'],      onupdate='CASCADE', ondelete='RESTRICT'),
        #ForeignKeyConstraint(['production_id'],     ['production.id'], onupdate='CASCADE', ondelete='CASCADE'),
        #UniqueConstraint('production_id', 'image_id', 'band', 'x', 'y'),
        )
    # Keys
    id = Column(                 BigInteger,   nullable=False) #Unique identifier"),
    #production_id = Column(       Integer,      nullable=False) #Production number"),
    image_id = Column(            BigInteger,   nullable=False) #CCD image number"),
    #Column('global_object_id = Column(    BigInteger,   nullable=True,  comment="Global Object number"),
    insert_date = Column(         DateTime,     nullable=False, default=func.current_timestamp()) #Timestamp of insertion"
    band = Column(                String(4),    nullable=False) #Band name"),
    # Fields
    background = Column(          Float(24),    nullable=False) #Background at centroid position [BACKGROUND] (count)"),
    class_star = Column(          Float(24),    nullable=False) #Star-Galaxy classifier output [CLASS_STAR]"),
    spread_model = Column(        Float(24),    nullable=True) #Spread parameter from model-fitting [SPREAD_MODEL]"),
    spreaderr_model = Column(     Float(24),    nullable=True) #Spread parameter error from model-fitting [SPREADERR_MODEL]"),
    flux_auto = Column(           Float(24),    nullable=False) #Flux within a Kron-like elliptical aperture [FLUX_AUTO] (count)"),
    flux_err_auto = Column(       Float(24),    nullable=False) #RMS error for AUTO flux  [FLUXERR_AUTO] (count)"),
    flux_psf = Column(            Float(24),    nullable=False) #Flux fitted using the PSF model (count)"),
    flux_err_psf = Column(        Float(24),    nullable=False) #RMS error for PSF flux"),
    flux_model = Column(          Float(24),    nullable=False) #Flux fitted using a galaxy model (count)"),
    flux_err_model = Column(      Float(24),    nullable=False) #RMS error for the model flux"),
    flags = Column(               SmallInteger, nullable=False) #Extraction flags [FLAGS]"),
    elongation = Column(          Float(24),    nullable=False) #A_IMAGE/B_IMAGE  [ELONGATION]"),
    dec = Column(                 Float(53),    nullable=False) #Windowed declination of barycenter (J2000) [DELTAWIN_J2000] (deg)"),
    ra = Column(                  Float(53),    nullable=False) #Windowed right ascension of barycenter (J2000) [ALPHAWIN_J2000] (deg)"),
    x = Column(                   Float(24),    nullable=False) #Windowed position estimate along x [XWIN_IMAGE] (pix)"),
    y = Column(                   Float(24),    nullable=False) #Windowed position estimate along y [YWIN_IMAGE] (pix)"),
    zp_offset = Column(           Float(24),    nullable=False) #Offset zeropoint magnitude for unexpected corrections"),
    #relationships
    #production  = relationship('Production',       back_populates='detections')
    image       = relationship('Image',            back_populates="detections")
    global_objects =  relationship('Global_object',    back_populates="detections", secondary=lambda:Global_object_detections.__table__)
    # Documentation



Index('ik_detlocation', Detection.ra, Detection.dec)


    # Global Object Table
class Global_object(Base):
    __tablename__ = 'global_object'
    __table_args__ = (
        # Constraints
        PrimaryKeyConstraint('id'),
        ForeignKeyConstraint(['production_id'], ['production.id'], onupdate='CASCADE', ondelete='CASCADE'),
        UniqueConstraint('production_id', 'ra', 'dec'),
            )   
    # Keys
    id = Column(                  BigInteger,   nullable=False) #Unique identifier"),
    production_id = Column(       Integer,      nullable=False) #Production number"),
    # Fields
    ra = Column(                  Float(53),    nullable=False) #Right Ascension"),
    dec = Column(                 Float(53),    nullable=False) #Declination"),
    sigma_offset_ra = Column(     Float(24),    nullable=False) #Sigma offset at RA"),
    sigma_offset_dec = Column(    Float(24),    nullable=False) #Sigma offset at DEC"),
        # Relationships
        
    production = relationship('Production',       back_populates="global_objects")
    detections = relationship('Detection',        back_populates="global_objects", secondary=lambda: Global_object_detections.__table__)
    coadd_objects = relationship('Coadd_object',     back_populates="global_object")
    
        # Documentation 
        #Links multiple detections to unique global objects.",

Index('ik_globallocation', Global_object.ra, Global_object.dec)

    # Global Object Detections Table
class Global_object_detections(Base):
    __tablename__ = 'global_object_detections'
    __table_args__ = (
        # Constraints
        PrimaryKeyConstraint('global_object_id', 'detection_id'),
        UniqueConstraint('detection_id', 'global_object_id'),
        ForeignKeyConstraint(['detection_id'],     ['detection.id'],     onupdate='CASCADE', ondelete='CASCADE'),
        ForeignKeyConstraint(['global_object_id'], ['global_object.id'], onupdate='CASCADE', ondelete='CASCADE'),
            )
        # Keys
    global_object_id =Column(   Integer,     nullable=False) #Global object identifier"),
    detection_id = Column(  BigInteger,  nullable=False) #Detection identifier"),
        # Constraints

        # Documentation
        #Links multiple detections to unique global objects.",
    # Coadd Object Table
class Coadd_object(Base):
    __tablename__ = 'coadd_object'
    __table_args__ = (
        # Constraints
        PrimaryKeyConstraint('id'),
        ForeignKeyConstraint(['production_id'],     ['production.id'], onupdate='CASCADE', ondelete='CASCADE'),
        ForeignKeyConstraint(['global_object_id'],  ['global_object.id'],      onupdate='CASCADE', ondelete='RESTRICT'),
        UniqueConstraint('production_id', 'ra', 'dec'),
            )    
        # Keys
    id = Column(                    BigInteger,   nullable=False) #Unique identifier"),
    production_id = Column(         Integer,      nullable=False) #Production number"),
    #template_image_tile_id = Column( BigInteger,   nullable=False) #CCD image number"),
    global_object_id = Column(      BigInteger,   nullable=True)
    insert_date = Column(           DateTime,     nullable=False) #Timestamp of insertion", default=func.current_timestamp()),
    # Fields
    ra = Column(                    Float(53),    nullable=False) #Windowed right ascension of barycenter (J2000) [ALPHAWIN_J2000] (deg)"),
    dec = Column(                   Float(53),    nullable=False) #Windowed declination of barycenter (J2000) [DELTAWIN_J2000] (deg)"),
    class_star = Column(            Float(24),    nullable=True) #Star-Galaxy classifier output [CLASS_STAR]"),
    spread_model = Column(          Float(24),    nullable=True) #Spread parameter from model-fitting [SPREAD_MODEL]"),
    spreaderr_model = Column(       Float(24),    nullable=True) #Spread parameter error from model-fitting [SPREADERR_MODEL]"),
    mag_model_u = Column(           Float(24),    nullable=True) #Magnitude fitted using a galaxy model"),
    mag_err_model_u = Column(       Float(24),    nullable=True) #RMS error for the model Magnitude"),
    mag_model_g = Column(           Float(24),    nullable=True) #Magnitude fitted using a galaxy model"),
    mag_err_model_g = Column(       Float(24),    nullable=True) #RMS error for the model Magnitude"),
    mag_model_r = Column(           Float(24),    nullable=True) #Magnitude fitted using a galaxy model"),
    mag_err_model_r = Column(       Float(24),    nullable=True) #RMS error for the model Magnitude"),
    mag_model_i = Column(           Float(24),    nullable=True) #Magnitude fitted using a galaxy model"),
    mag_err_model_i = Column(       Float(24),    nullable=True) #RMS error for the model Magnitude"),
    mag_model_z = Column(           Float(24),    nullable=True) #Magnitude fitted using a galaxy model"),
    mag_err_model_z = Column(       Float(24),    nullable=True) #RMS error for the model Magnitude"),
    mag_model_Y = Column(           Float(24),    nullable=True) #Magnitude fitted using a galaxy model"),
    mag_err_model_Y = Column(       Float(24),    nullable=True) #RMS error for the model Magnitude"),
    mag_model_n01 = Column(         Float(24),    nullable=True) #Magnitude fitted using a galaxy model"),
    mag_err_model_n01 = Column(     Float(24),    nullable=True) #RMS error for the model Magnitude"),
    mag_model_n02 = Column(         Float(24),    nullable=True) #Magnitude fitted using a galaxy model"),
    mag_err_model_n02 = Column(     Float(24),    nullable=True) #RMS error for the model Magnitude"),
    mag_model_n03 = Column(         Float(24),    nullable=True) #Magnitude fitted using a galaxy model"),
    mag_err_model_n03 = Column(     Float(24),    nullable=True) #RMS error for the model Magnitude"),
    mag_model_n04 = Column(         Float(24),    nullable=True) #Magnitude fitted using a galaxy model"),
    mag_err_model_n04 = Column(     Float(24),    nullable=True) #RMS error for the model Magnitude"),
    mag_model_n05 = Column(         Float(24),    nullable=True) #Magnitude fitted using a galaxy model"),
    mag_err_model_n05 = Column(     Float(24),    nullable=True) #RMS error for the model Magnitude"),
    mag_model_n06 = Column(         Float(24),    nullable=True) #Magnitude fitted using a galaxy model"),
    mag_err_model_n06 = Column(     Float(24),    nullable=True) #RMS error for the model Magnitude"),
    mag_model_n07 = Column(         Float(24),    nullable=True) #Magnitude fitted using a galaxy model"),
    mag_err_model_n07 = Column(     Float(24),    nullable=True) #RMS error for the model Magnitude"),
    mag_model_n08 = Column(         Float(24),    nullable=True) #Magnitude fitted using a galaxy model"),
    mag_err_model_n08 = Column(     Float(24),    nullable=True) #RMS error for the model Magnitude"),
    mag_model_n09 = Column(         Float(24),    nullable=True) #Magnitude fitted using a galaxy model"),
    mag_err_model_n09 = Column(     Float(24),    nullable=True) #RMS error for the model Magnitude"),
    mag_model_n10 = Column(         Float(24),    nullable=True) #Magnitude fitted using a galaxy model"),
    mag_err_model_n10 = Column(     Float(24),    nullable=True) #RMS error for the model Magnitude"),
    mag_model_n11 = Column(         Float(24),    nullable=True) #Magnitude fitted using a galaxy model"),
    mag_err_model_n11 = Column(     Float(24),    nullable=True) #RMS error for the model Magnitude"),
    mag_model_n12 = Column(         Float(24),    nullable=True) #Magnitude fitted using a galaxy model"),
    mag_err_model_n12 = Column(     Float(24),    nullable=True) #RMS error for the model Magnitude"),
    mag_model_n13 = Column(         Float(24),    nullable=True) #Magnitude fitted using a galaxy model"),
    mag_err_model_n13 = Column(     Float(24),    nullable=True) #RMS error for the model Magnitude"),
    mag_model_n14 = Column(         Float(24),    nullable=True) #Magnitude fitted using a galaxy model"),
    mag_err_model_n14 = Column(     Float(24),    nullable=True) #RMS error for the model Magnitude"),
    mag_model_n15 = Column(         Float(24),    nullable=True) #Magnitude fitted using a galaxy model"),
    mag_err_model_n15 = Column(     Float(24),    nullable=True) #RMS error for the model Magnitude"),
    mag_model_n16 = Column(         Float(24),    nullable=True) #Magnitude fitted using a galaxy model"),
    mag_err_model_n16 = Column(     Float(24),    nullable=True) #RMS error for the model Magnitude"),
    mag_model_n17 = Column(         Float(24),    nullable=True) #Magnitude fitted using a galaxy model"),
    mag_err_model_n17 = Column(     Float(24),    nullable=True) #RMS error for the model Magnitude"),
    mag_model_n18 = Column(         Float(24),    nullable=True) #Magnitude fitted using a galaxy model"),
    mag_err_model_n18 = Column(     Float(24),    nullable=True) #RMS error for the model Magnitude"),
    mag_model_n19 = Column(         Float(24),    nullable=True) #Magnitude fitted using a galaxy model"),
    mag_err_model_n19 = Column(     Float(24),    nullable=True) #RMS error for the model Magnitude"),
    mag_model_n20 = Column(         Float(24),    nullable=True) #Magnitude fitted using a galaxy model"),
    mag_err_model_n20 = Column(     Float(24),    nullable=True) #RMS error for the model Magnitude"),
    mag_model_n21 = Column(         Float(24),    nullable=True) #Magnitude fitted using a galaxy model"),
    mag_err_model_n21 = Column(     Float(24),    nullable=True) #RMS error for the model Magnitude"),
    mag_model_n22 = Column(         Float(24),    nullable=True) #Magnitude fitted using a galaxy model"),
    mag_err_model_n22 = Column(     Float(24),    nullable=True) #RMS error for the model Magnitude"),
    mag_model_n23 = Column(         Float(24),    nullable=True) #Magnitude fitted using a galaxy model"),
    mag_err_model_n23 = Column(     Float(24),    nullable=True) #RMS error for the model Magnitude"),
    mag_model_n24 = Column(         Float(24),    nullable=True) #Magnitude fitted using a galaxy model"),
    mag_err_model_n24 = Column(     Float(24),    nullable=True) #RMS error for the model Magnitude"),
    mag_model_n25 = Column(         Float(24),    nullable=True) #Magnitude fitted using a galaxy model"),
    mag_err_model_n25 = Column(     Float(24),    nullable=True) #RMS error for the model Magnitude"),
    mag_model_n26 = Column(         Float(24),    nullable=True) #Magnitude fitted using a galaxy model"),
    mag_err_model_n26 = Column(     Float(24),    nullable=True) #RMS error for the model Magnitude"),
    mag_model_n27 = Column(         Float(24),    nullable=True) #Magnitude fitted using a galaxy model"),
    mag_err_model_n27 = Column(     Float(24),    nullable=True) #RMS error for the model Magnitude"),
    mag_model_n28 = Column(         Float(24),    nullable=True) #Magnitude fitted using a galaxy model"),
    mag_err_model_n28 = Column(     Float(24),    nullable=True) #RMS error for the model Magnitude"),
    mag_model_n29 = Column(         Float(24),    nullable=True) #Magnitude fitted using a galaxy model"),
    mag_err_model_n29 = Column(     Float(24),    nullable=True) #RMS error for the model Magnitude"),
    mag_model_n30 = Column(         Float(24),    nullable=True) #Magnitude fitted using a galaxy model"),
    mag_err_model_n30 = Column(     Float(24),    nullable=True) #RMS error for the model Magnitude"),
    mag_model_n31 = Column(         Float(24),    nullable=True) #Magnitude fitted using a galaxy model"),
    mag_err_model_n31 = Column(     Float(24),    nullable=True) #RMS error for the model Magnitude"),
    mag_model_n32 = Column(         Float(24),    nullable=True) #Magnitude fitted using a galaxy model"),
    mag_err_model_n32 = Column(     Float(24),    nullable=True) #RMS error for the model Magnitude"),
    mag_model_n33 = Column(         Float(24),    nullable=True) #Magnitude fitted using a galaxy model"),
    mag_err_model_n33 = Column(     Float(24),    nullable=True) #RMS error for the model Magnitude"),
    mag_model_n34 = Column(         Float(24),    nullable=True) #Magnitude fitted using a galaxy model"),
    mag_err_model_n34 = Column(     Float(24),    nullable=True) #RMS error for the model Magnitude"),
    mag_model_n35 = Column(         Float(24),    nullable=True) #Magnitude fitted using a galaxy model"),
    mag_err_model_n35 = Column(     Float(24),    nullable=True) #RMS error for the model Magnitude"),
    mag_model_n36 = Column(         Float(24),    nullable=True) #Magnitude fitted using a galaxy model"),
    mag_err_model_n36 = Column(     Float(24),    nullable=True) #RMS error for the model Magnitude"),
    mag_model_n37 = Column(         Float(24),    nullable=True) #Magnitude fitted using a galaxy model"),
    mag_err_model_n37 = Column(     Float(24),    nullable=True) #RMS error for the model Magnitude"),
    mag_model_n38 = Column(         Float(24),    nullable=True) #Magnitude fitted using a galaxy model"),
    mag_err_model_n38 = Column(     Float(24),    nullable=True) #RMS error for the model Magnitude"),
    mag_model_n39 = Column(         Float(24),    nullable=True) #Magnitude fitted using a galaxy model"),
    mag_err_model_n39 = Column(     Float(24),    nullable=True) #RMS error for the model Magnitude"),
    mag_model_n40 = Column(         Float(24),    nullable=True) #Magnitude fitted using a galaxy model"),
    mag_err_model_n40 = Column(     Float(24),    nullable=True) #RMS error for the model Magnitude"),
    mag_model_n41 = Column(         Float(24),    nullable=True) #Magnitude fitted using a galaxy model"),
    mag_err_model_n41 = Column(     Float(24),    nullable=True) #RMS error for the model Magnitude"),
    mag_model_n42 = Column(         Float(24),    nullable=True) #Magnitude fitted using a galaxy model"),
    mag_err_model_n42 = Column(     Float(24),    nullable=True) #RMS error for the model Magnitude"),
    flags = Column(                 SmallInteger, nullable=True) #Extraction flags [FLAGS]"),
    elongation = Column(            Float(24),    nullable=False) #A_IMAGE/B_IMAGE  [ELONGATION]"),
    
    #Reletionships
    production    = relationship('Production',       back_populates="coadd_objects")
    global_object = relationship('Global_object',    back_populates="coadd_objects")
    photo_zs      = relationship('Photo_z'     ,     back_populates="coadd_object")
    
# Documentation
#comment="Contains unique coadd objects extracted from coadd image tiles.",

Index('ik_coaddlocation', Coadd_object.production_id, Coadd_object.ra, Coadd_object.dec)

class Photo_z(Base):
    __tablename__ = 'photo_z'
    __table_args__ = (
        # Constraints
        PrimaryKeyConstraint('id'),
        ForeignKeyConstraint(['production_id'], ['production.id'], onupdate='CASCADE', ondelete='CASCADE'),
        ForeignKeyConstraint(['coadd_object_id'], ['coadd_object.id'], onupdate='CASCADE', ondelete='CASCADE'),
        UniqueConstraint('production_id', 'coadd_object_id', 'code'),
            )     
    # Keys
    id = Column(                    BigInteger,   nullable=False) #Unique identifier"),
    production_id = Column(         Integer,      nullable=False) #Production number"),
    coadd_object_id = Column(       BigInteger,   nullable=False) #Production number"),
    # Fields
    z = Column(                     Float(24),    nullable=False) #Photometric Redshift"),
    z_err = Column(                 Float(24),    nullable=False) #Photometric Redshift Error"),
    odds = Column(                  Float(24),    nullable=False) #Redshift odds"),
    pdf = Column(                   String(2000), nullable=False) #Redshift Probability Density Function"),
    sed_type = Column(              Integer,      nullable=False) #Spectral Energy Distribution type"),
    code = Column(                  Integer,      nullable=False) #Code used to determine redshift"),
    # Relationships

    production    = relationship('Production',       back_populates="photo_zs")
    coadd_object = relationship('Coadd_object',     back_populates="photo_zs")

    # Documentation
    #comment="Links multiple detections to unique global objects.",

Index('ik_photozredshift', Photo_z.production_id, Photo_z.z)
    

class SDSS(Base):
    __tablename__ = 'sdss'
    __table_args__ = (
        # Constraints
        PrimaryKeyConstraint('id'),
            )     
        
    # Keys
    id = Column(        BigInteger,   nullable=False) #Unique identifier. Same as Unique SDSS identifier composed from [skyVersion,rerun,run,camcol,field,obj]"),
    # Fields
    ra = Column(        Float(53),    nullable=False) #Right Ascension of the object [RAJ2000] (deg)"),
    dec = Column(       Float(53),    nullable=False) #Declination of the object [DEJ2000] (deg)"),
    ra_err = Column(    Float(24),    nullable=False) #Error in RA (arcsec)"),
    dec_err = Column(   Float(24),    nullable=False) #Error in DEC (arcsec)"),
    mag_u = Column(     Float(24),    nullable=False) #Model magnitude in u filter"),
    mag_g = Column(     Float(24),    nullable=False) #Model magnitude in g filter"),
    mag_r = Column(     Float(24),    nullable=False) #Model magnitude in r filter"),
    mag_i = Column(     Float(24),    nullable=False) #Model magnitude in i filter"),
    mag_z = Column(     Float(24),    nullable=False) #Model magnitude in z filter"),
    mag_err_u = Column( Float(24),    nullable=False) #Mean error on umag"),
    mag_err_g = Column( Float(24),    nullable=False) #Mean error on gmag"),
    mag_err_r = Column( Float(24),    nullable=False) #Mean error on rmag"),
    mag_err_i = Column( Float(24),    nullable=False) #Mean error on imag"),
    mag_err_z = Column( Float(24),    nullable=False) #Mean error on zmag"),
    type = Column(      SmallInteger, nullable=False) #Morphological type classification of the object. 0 UNKOWN, 1 COSMIC_RAY, 2 DEFECT, 3 GALAXY, 4 GHOST, 5 KNOWNOBJ, 6 STAR, 7 TRAIL, 8 SKY, 9 NOTATYPE"),
    prob_psf = Column(  Float(24),    nullable=False) #Probability that the object is a star. Currently 0 if type == 3 (galaxy), 1 if type == 6 (star)."),
    u_s = Column(       Float(24),    nullable=False) #[0,1] 0=notStar, 1=Star in u band (probPSF_u)"),
    g_s = Column(       Float(24),    nullable=False) #[0,1] 0=notStar, 1=Star in g band (probPSF_g)"),
    r_s = Column(       Float(24),    nullable=False) #[0,1] 0=notStar, 1=Star in r band (probPSF_r)"),
    i_s = Column(       Float(24),    nullable=False) #[0,1] 0=notStar, 1=Star in i band (probPSF_i)"),
    z_s = Column(       Float(24),    nullable=False) #[0,1] 0=notStar, 1=Star in z band (probPSF_z)"),
    # Constraints
    PrimaryKeyConstraint('id'),
    # Documentation
    #comment="SDSS",

Index('ik_sdsslocation', SDSS.ra, SDSS.dec)
    
    # USNO External Table
class USNO(Base):
    __tablename__ = 'usno'
    __table_args__ = (
        # Constraints
        PrimaryKeyConstraint('id'),
            )   
        # Keys
    id = Column(        BigInteger,   nullable=False) #Unique identifier"),
    # Fields
    ra = Column(        Float(53),    nullable=False) #Right Ascension of the object (deg)"),
    dec = Column(       Float(53),    nullable=False) #Declination of the object (deg)"),
    ra_err = Column(    Float(24),    nullable=False) #Error in RA (arcsec)"),
    dec_err = Column(   Float(24),    nullable=False) #Error in DEC (arcsec)"),
    mag_r = Column(     Float(24),    nullable=False) #Model magnitude in r filter"),
    mag_b = Column(     Float(24),    nullable=False) #Model magnitude in b filter"),
    mag_err_r = Column( Float(24),    nullable=False) #Mean error on mag_r"),
    mag_err_b = Column( Float(24),    nullable=False) #Mean error on mag_b"),

Index('ik_usnolocation', USNO.ra, USNO.dec)
  
    
    ###########    EXTERNAL Tables    ###########
    
    # SDSS External Table

    
    ###########    SIMULATION Tables    ###########
    
    # Truth Objects table
class Truth_Object(Base):
    __tablename__ = 'truth_object'
    __table_args__ = (
        # Constraints
        PrimaryKeyConstraint('id'),
        ForeignKeyConstraint(['production_id'], ['production.id'], onupdate='CASCADE', ondelete='CASCADE'),
            ) 
    # Keys
    id = Column(               BigInteger,   nullable=False) #Unique identifier"),
    # Fields
    production_id = Column(    Integer,      nullable=False) #Production number"),
    stargalaxy = Column(       Boolean,      nullable=False) #0(False):galaxy, 1(True):star"),
    sed_type = Column(         Float(24),    nullable=False) #SED type"),
    sed_em_line = Column(      Float(24),    nullable=False) #SED emission line type"),
    ra = Column(               Float(53),    nullable=False) #Truth Right AscensSion position"),
    dec = Column(              Float(53),    nullable=False) #Truth Declination position"),
    mag_u = Column(            Float(24),    nullable=False) #Truth simulated magnitude for filter u"),
    mag_g = Column(            Float(24),    nullable=False) #Truth simulated magnitude for filter g"),
    mag_r = Column(            Float(24),    nullable=False) #Truth simulated magnitude for filter r"),
    mag_i = Column(            Float(24),    nullable=False) #Truth simulated magnitude for filter i"),
    mag_z = Column(            Float(24),    nullable=False) #Truth simulated magnitude for filter z"),
    mag_y = Column(            Float(24),    nullable=False) #Truth simulated magnitude for filter y"),
    mag_NB455 = Column(          Float(24),    nullable=False) #Truth simulated magnitude for PAU narrow filter 01"),
    mag_NB465 = Column(          Float(24),    nullable=False) #Truth simulated magnitude for PAU narrow filter 02"),
    mag_NB475 = Column(          Float(24),    nullable=False) #Truth simulated magnitude for PAU narrow filter 03"),
    mag_NB485 = Column(          Float(24),    nullable=False) #Truth simulated magnitude for PAU narrow filter 04"),
    mag_NB495 = Column(          Float(24),    nullable=False) #Truth simulated magnitude for PAU narrow filter 05"),
    mag_NB505 = Column(          Float(24),    nullable=False) #Truth simulated magnitude for PAU narrow filter 06"),
    mag_NB515 = Column(          Float(24),    nullable=False) #Truth simulated magnitude for PAU narrow filter 07"),
    mag_NB525 = Column(          Float(24),    nullable=False) #Truth simulated magnitude for PAU narrow filter 08"),
    mag_NB535 = Column(          Float(24),    nullable=False) #Truth simulated magnitude for PAU narrow filter 09"),
    mag_NB545 = Column(          Float(24),    nullable=False) #Truth simulated magnitude for PAU narrow filter 10"),
    mag_NB555 = Column(          Float(24),    nullable=False) #Truth simulated magnitude for PAU narrow filter 11"),
    mag_NB565 = Column(          Float(24),    nullable=False) #Truth simulated magnitude for PAU narrow filter 12"),
    mag_NB575 = Column(          Float(24),    nullable=False) #Truth simulated magnitude for PAU narrow filter 13"),
    mag_NB585 = Column(          Float(24),    nullable=False) #Truth simulated magnitude for PAU narrow filter 14"),
    mag_NB595 = Column(          Float(24),    nullable=False) #Truth simulated magnitude for PAU narrow filter 15"),
    mag_NB605 = Column(          Float(24),    nullable=False) #Truth simulated magnitude for PAU narrow filter 16"),
    mag_NB615 = Column(          Float(24),    nullable=False) #Truth simulated magnitude for PAU narrow filter 17"),
    mag_NB625 = Column(          Float(24),    nullable=False) #Truth simulated magnitude for PAU narrow filter 18"),
    mag_NB635 = Column(          Float(24),    nullable=False) #Truth simulated magnitude for PAU narrow filter 19"),
    mag_NB645 = Column(          Float(24),    nullable=False) #Truth simulated magnitude for PAU narrow filter 20"),
    mag_NB655 = Column(          Float(24),    nullable=False) #Truth simulated magnitude for PAU narrow filter 21"),
    mag_NB665 = Column(          Float(24),    nullable=False) #Truth simulated magnitude for PAU narrow filter 22"),
    mag_NB675 = Column(          Float(24),    nullable=False) #Truth simulated magnitude for PAU narrow filter 23"),
    mag_NB685 = Column(          Float(24),    nullable=False) #Truth simulated magnitude for PAU narrow filter 24"),
    mag_NB695 = Column(          Float(24),    nullable=False) #Truth simulated magnitude for PAU narrow filter 25"),
    mag_NB705 = Column(          Float(24),    nullable=False) #Truth simulated magnitude for PAU narrow filter 26"),
    mag_NB715 = Column(          Float(24),    nullable=False) #Truth simulated magnitude for PAU narrow filter 27"),
    mag_NB725 = Column(          Float(24),    nullable=False) #Truth simulated magnitude for PAU narrow filter 28"),
    mag_NB735 = Column(          Float(24),    nullable=False) #Truth simulated magnitude for PAU narrow filter 29"),
    mag_NB745 = Column(          Float(24),    nullable=False) #Truth simulated magnitude for PAU narrow filter 30"),
    mag_NB755 = Column(          Float(24),    nullable=False) #Truth simulated magnitude for PAU narrow filter 31"),
    mag_NB765 = Column(          Float(24),    nullable=False) #Truth simulated magnitude for PAU narrow filter 32"),
    mag_NB775 = Column(          Float(24),    nullable=False) #Truth simulated magnitude for PAU narrow filter 33"),
    mag_NB785 = Column(          Float(24),    nullable=False) #Truth simulated magnitude for PAU narrow filter 34"),
    mag_NB795 = Column(          Float(24),    nullable=False) #Truth simulated magnitude for PAU narrow filter 35"),
    mag_NB805 = Column(          Float(24),    nullable=False) #Truth simulated magnitude for PAU narrow filter 36"),
    mag_NB815 = Column(          Float(24),    nullable=False) #Truth simulated magnitude for PAU narrow filter 37"),
    mag_NB825 = Column(          Float(24),    nullable=False) #Truth simulated magnitude for PAU narrow filter 38"),
    mag_NB835 = Column(          Float(24),    nullable=False) #Truth simulated magnitude for PAU narrow filter 39"),
    mag_NB845 = Column(          Float(24),    nullable=False) #Truth simulated magnitude for PAU narrow filter 40"),
    mag_n41 = Column(          Float(24),    nullable=True) #Truth simulated magnitude for PAU narrow filter 41"),
    mag_n42 = Column(          Float(24),    nullable=True) #Truth simulated magnitude for PAU narrow filter 42"),
    bulge2flux_ratio = Column( Float(24),    nullable=True) #Galaxy only - Bulge to total flux ratio"),
    bulge_length = Column(     Float(24),    nullable=True) #Galaxy only - Bulge length scale (arcsec)"),
    bulge_axis_ratio = Column( Float(24),    nullable=True) #Galaxy only - Bulge projected axis ratio"),
    bulge_angle = Column(      Float(24),    nullable=True) #Galaxy only - Bulge position angle"),
    disk_length = Column(      Float(24),    nullable=True) #Galaxy only - Disk length scale (arcsec)"),
    disk_axis_ratio = Column(  Float(24),    nullable=True) #Galaxy only - Disk projected axis ratio"),
    disk_angle = Column(       Float(24),    nullable=True) #Galaxy only - Disk position angle"),
    redshift = Column(         Float(24),    nullable=True) #Galaxy only - Truth redshift"),
    # Constraints
    production    = relationship('Production',       back_populates="truth_objects")
    # Documentation
    #comment="The truth object table contains all objects used in the simulated images for all PAU data challenges.\n"
    #"Each row corresponds to a truth object and for this reason all parameters do not have associated errors.",
Index('ik_truthlocation', Truth_Object.ra, Truth_Object.dec)
    
    # Target
class Target(Base):
    __tablename__ = 'target'
    __table_args__ = (
        # Constraints
        PrimaryKeyConstraint('id'),
        ForeignKeyConstraint(['production_id'], ['production.id'], onupdate='CASCADE', ondelete='CASCADE'),
            ) 
    # Keys
    id = Column(                Integer,      nullable=False) #Unique identifier"),
    production_id = Column(     Integer,      nullable=False) #Production number"),
    # Fields
    exp_num = Column(           Integer,      nullable=True) #  comment="Camera exposure unique identifier"),
    field = Column(             String(16),   nullable=True) #  comment="Field Name"),
    ra = Column(                Float(53),    nullable=True) #Target RA"),
    dec = Column(               Float(53),    nullable=True) #Target DEC"),
    filtertray = Column(        String(16),   nullable=True) #Filter Tray name"),
    kind = Column(              Enum('BIAS', 'FLAT', 'TARGET', name='target_kind'), nullable=False) #Target kind"),
    status = Column(            Enum('PLANNED', 'SCHEDULED', 'SIMULATED', 'EXPOSED', name='target_status'),nullable=False) #  comment="Target status"),
    rjd_obs = Column(           Float(53),    nullable=True) #Observation Reduced Modified Julian Day"),
    exptime = Column(           Float(24),    nullable=True) # comment="Exposure Time",
    airmass = Column(           Float(24),    nullable=True) # comment="Airmass",
    seeing  = Column(           Float(24),    nullable=True) # comment="Seeing",
    moon_mag = Column(          Float(24),    nullable=True) # comment="Moon Magnitude",
    moon_phase = Column(        Float(24),    nullable=True) # comment="Moon Phase",
    moon_distance =  Column(    Float(24),    nullable=True) # comment="Moon Distance (degrees)",
    moon_set = Column(          Boolean,      nullable=True) # comment="Moon is set?",
    ecliptic_dist = Column(     Float(24),    nullable=True) # comment="Distance to ecliptic (degrees)",
    ecliptic_zodiacal = Column( Float(24),  nullable=True) # comment="Ecliptic Zodiacal Light",
    ecliptic_airglow = Column(  Float(24),  nullable=True) # comment="Ecliptic Airglow Light",
    
    # Relationships
    bkg_mags   = relationship('Bkg_mag',    back_populates="target")
    production = relationship('Production', back_populates="targets")
    # Documentation
    #comment="target",

Index('ik_targetlocation', Target.ra, Target.dec)

class Bkg_mag(Base):
    __tablename__ = 'bkg_mag'
    __table_args__ = (
        # Constraints
        PrimaryKeyConstraint('target_id', 'filter'),
        ForeignKeyConstraint(['target_id'], ['target.id'], ondelete='CASCADE')
        )
    # Keys
    target_id = Column( Integer,      nullable=False) # comment="Unique identifier",
    filter    = Column( String(16),   nullable=False) # comment="Filter Name"),
    # Fields
    mag       = Column( Float(24),    nullable=False) # comment="Magnitude value"),
    # Documentation
    #comment="Background Magnitude",
    
    # Relationships
    target = relationship('Target', back_populates="bkg_mags")
    
    ###########    OPERATION Tables    ###########
    
    
    # Job Table
# class Job_pau(brownthrower.model.Job) :
#         
#         quality_controls = relationship('Quality_control',  back_populates="job")   
#     
#     # Documentation
#     #comment="Job",
#     
#     # Quality Control Table
class Quality_control(Base):
    __tablename__ = 'quality_control'
    __table_args__ = (
        # Constraints
        PrimaryKeyConstraint('id'),
        
            )
    # Keys
    id = Column(                Integer,      nullable=False) #Unique identifier"),
    # Fields
    job_id = Column(            Integer,      nullable=False) #Job identifier"),
    ref = Column(               String(24),   nullable=False) #QC Reference code"),
    check_name = Column(        String(64),   nullable=False) #QC check name"),
    min_value = Column(         Float(24),    nullable=False) #Minimum Range value"),
    max_value = Column(         Float(24),    nullable=False) #Maximum Range value"),
    value = Column(             Float(24),    nullable=False) #Measured value"),
    units = Column(             String(24),   nullable=False) #Measurement units"),
    qc_pass = Column(           Boolean,      nullable=False) #QC pass?"),
    time = Column(              DateTime,     nullable=False,default=func.current_timestamp() ) #comment="Timestamp of creation", ),
    plot_file = Column(         Text,         nullable=True) #  comment="Plot file full name and path"),
    
    # Relationships
#     job           = relationship('Job_pau',              back_populates="quality_controls")
    
    #comment="Quality Control",

    
@compiles(BigInteger, "sqlite")
def compile_biginteger_sqlite(type_, compiler, **kw):
    """
    SQLite does not allow autoincremented BIGINT primary keys.
    This method rewrites those colums as INTEGER to overcome this limitation.
    See: http://www.sqlite.org/faq.html#q1
    """
    return "INTEGER"

def add_criteria(query, criteria):
    q = query
    
    for table_name, columns in criteria.iteritems():
        for column_name, filtr in columns.iteritems():
            if isinstance(filtr, list):
                q = q.filter(tables[table_name].columns[column_name].in_(filtr))
            elif isinstance(filtr, dict):
                for operator, value in filtr.iteritems():
                    q = q.filter(tables[table_name].columns[column_name].op(operator)(value))
            else:
                q = q.filter(tables[table_name].columns[column_name] == filtr)
    return q



def opt_parser(argv):
    parser = optparse.OptionParser(usage="Usage: %prog -u URL")
    parser.add_option("-u", "--url", dest="url",
                      default="sqlite:///",
                      help="URL for database connection [default: %default].")
    parser.add_option("-f", "--force", dest="force",
                      action="store_true", default=False,
                      help="Recreate database without confirmation [default: %default].")
    
    (options, args) = parser.parse_args(argv)
    if not options.url:
        parser.error("You must specify a database URL (-u).")
    
    return (options, args)

def main(argv = None):
    if not argv:
        argv = sys.argv
    
    (options, ) = opt_parser(argv)[:1]

    init(options.url)
    
    if not options.force:
        answer = raw_input("The database will be ERASED and recreated. Do you want to proceed? (y/N): ")
        if answer.upper() != "Y":
            return
    
    recreate()

def recreate():
    # Drop and recreate the database
    metadata.drop_all()
    #metadata.create_all()
    if metadata.bind.url.database == 'paudb' and metadata.bind.url.drivername == 'postgresql':
        metadata.create_comments()
        metadata.set_permissions()

if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
