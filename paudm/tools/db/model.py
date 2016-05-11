#!/usr/bin/env python
# -*- coding: utf-8 -*-
import datetime
import hashlib
import optparse
import yaml
 
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
    

    twophase = False
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
   
    # Comment:
    # Contains the list of Observation Sets registered in the database.

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
    pipeline            = Column(String(32), nullable=True)   # Pipeline Name [pixelsim, nightly, memba, analysis]
    release             = Column(String(64), nullable=False)  # Major release name [TESTX, DRX]
    software_version    = Column(String(32), nullable=False)  # Package version [DCX, vX.X]
    _comments           = Column('comments',  Text,       nullable=True)
    job_id              = Column(Integer,    nullable=False)  # id of the job that generates the current production

    #Relationships
    mosaics        = relationship('Mosaic',         back_populates="production")
    global_objects = relationship('Global_object',  back_populates="production")
    coadd_objects  = relationship('Coadd_object',   back_populates="production")
    forced_aperture_coadds = relationship('ForcedApertureCoadd',  back_populates="production")
    forced_apertures = relationship('ForcedAperture',  back_populates="production")
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
    
    @hybrid_property
    def comments(self):
        return self._comments
    
    def get_comments(self):
        c = self.comments
        if c is None:
            return []
        else:
            return yaml.safe_load(c)
    
    def add_comment(self, value):
        c = self.get_comments()
        
        tnow = datetime.datetime.strftime(datetime.datetime.now(), '%Y-%m-%d %H:%M:%S')
        c.append({
            'timestamp' : tnow,
            'comment' : value
        })
        self._comments = yaml.safe_dump(c, default_flow_style=False)

    # Comment:
    # Tracks the different processing production runs for all pipelines.


class Mosaic(Base): 
    """
    Defines a standard Mosaic object.
    """

    __tablename__ = 'mosaic'
    __table_args__ = (
        # Primary key
        PrimaryKeyConstraint('id'),

        # Foreign key
        ForeignKeyConstraint(['production_id'], ['production.id'], onupdate='CASCADE', ondelete='CASCADE'),
        ForeignKeyConstraint(['obs_set_id'],        ['obs_set.id'], onupdate='CASCADE', ondelete='CASCADE'),
        UniqueConstraint('archivepath', 'filename'),
        UniqueConstraint('production_id', 'obs_set_id', 'kind', 'exp_num'),

    )
    # Columns
    id = Column(BigInteger, nullable=False)  # Identifier
    production_id = Column(Integer, nullable=False)  # Production identifier
    obs_set_id = Column(Integer, nullable=False)  # obs_set number
    filename = Column(String(128), nullable=False)  # File name
    comment = Column(Text, nullable=True)
    archivepath = Column(String(128), nullable=False)  # Path in the archive
    kind = Column(Enum('ARC', 'BIAS', 'DARK', 'FLASH', 'FLAT', 'FOCUS', 'GLANCE', 'PUPIL', 'SCIENCE',
                       'SCRATCH', 'SKY', 'TARGET','MBIAS','MFLAT','RED_SCI','RED_WEIGHT','RED_MASK',
                       name='mosaic_kind'), nullable=False)  # Mosaic image type
    exp_num = Column(Integer, nullable=False)  # Exposure number
    obs_title = Column(String(128), nullable=False)  # Observation title
    ra = Column(Float(53), nullable=False)  # Telescope Right Ascension pointing (deg)
    dec = Column(Float(53), nullable=False)  # Telescope Declination pointing (deg)
    equinox = Column(String(32), nullable=False)  # Equinox of telescope coordinates
    date_obs = Column(Date, nullable=False)  # Observation date
    time_obs = Column(Time, nullable=False)  # Observation time (Universal time)
    rjd_obs = Column(Float(53), nullable=False)  # Observation Reduced Modified Julian Day
    date_creat = Column(Date, nullable=False)  # File creation date
    time_creat = Column(Time, nullable=False)  # File creation time
    exp_time = Column(Float(24), nullable=False)  # Exposure time (s)
    airmass = Column(Float(24), nullable=True)  # airmass
    telfocus = Column(Float(24), nullable=True)  # Telescope Focus
    instrument = Column(String(32), nullable=False)  # Instrument name
    filtertray = Column(String(16), nullable=True)  # Filter Tray name
    filtertray_tmp = Column(Float(24), nullable=True)  # Filter Tray temperature (deg C)
    nextend = Column(SmallInteger, nullable=False)  # Filter Tray name
    guide_enabled = Column(Boolean, nullable=True)  # Number of extensions
    guide_fwhm = Column(Float(24), nullable=True)  # Seeing FWHM measured at guiding star (arcsec)
    guide_var = Column(Float(24), nullable=True)  # Flux variance measured at guiding star (counts)
    extinction = Column(Float(24), nullable=True)  # Extinction value measured in reduction for photometric calibration. Available in reduced image only.
    extinction_err = Column(Float(24), nullable=True)  # Extinction error value. Available in reduced image only.
    wind_spd = Column(Float(24), nullable=True)  # Wind speed (kph)
    wind_dir = Column(Float(24), nullable=True)  # Wind direction (deg)
    amb_temp = Column(Float(24), nullable=True)  # Ambient temperature (deg C)
    humidity = Column(Float(24), nullable=True)  # Ambient relative humidity (percent)
    pressure = Column(Float(24), nullable=True)  # Barometic pressure (mbar)
    eqa_1 = Column(Float(24), nullable=True)  # Quality Analysis check at observatory 1
    eqa_2 = Column(Float(24), nullable=True)  # Quality Analysis check at observatory 2
    eqa_3 = Column(Float(24), nullable=True)  # Quality Analysis check at observatory 3
    eqa_4 = Column(Float(24), nullable=True)  # Quality Analysis check at observatory 4
    eqa_5 = Column(Float(24), nullable=True)  # Quality Analysis check at observatory 5
    merged_mosaics = Column(Integer, nullable=True)  # Number of mosaics merged to form the actual mosaic (for masters)
    mean_psf_fwhm = Column(Float(24), nullable=True)  # Mean PSF FWHM measured. Available in reduced image only,

    # Relationships
    production = relationship('Production', back_populates="mosaics")
    obs_set = relationship('Obs_set', back_populates="mosaics")
    images = relationship('Image', back_populates="mosaic")

    @classmethod
    def criteria_query(cls, criteria={}):
        """
        run a query to retrieve Mosaic instances to be prestaged
        """
        q = session.query(cls).join(cls.production).join(cls.obs_set)
        q = add_criteria(q, criteria)
        return q

Index('ik_location', Mosaic.production_id, Mosaic.ra, Mosaic.dec)
Index('ik_rjdobs', Mosaic.production_id, Mosaic.rjd_obs)

    # Comment:
    # Contains the list of mosaic exposure images (raw and reduced).


class Image(Base):
    """
    Defines a standard Image object.

    """
    __tablename__ = 'image'
    __table_args__ = (
        PrimaryKeyConstraint('id'),
        ForeignKeyConstraint(['mosaic_id'], ['mosaic.id'], onupdate='CASCADE', ondelete='CASCADE'),
        UniqueConstraint('mosaic_id', 'ccd_num', 'amp_num'),
        UniqueConstraint('mosaic_id', 'image_num'),
        )
    # Keys
    id = Column(BigInteger, nullable=False)
    mosaic_id = Column(BigInteger, nullable=False)
    # Fields
    image_num = Column(SmallInteger, nullable=False)  # Extension number
    ccd_num = Column(SmallInteger, nullable=False)  # CCD Number"
    amp_num = Column(SmallInteger, nullable=False)  # Amplifier number (-1 for full CCD)
    filter = Column(String(8), nullable=True)   # Filter name"),
    wavelength = Column(Float(24), nullable=False)  # Wavelength at filter center (nm)
    waveband = Column(Float(24), nullable=False)  # Filter's Waveband
    gain = Column(Float(24), nullable=False)  # Detector gain at amplifier (e-/ADU)
    rdnoise = Column(Float(24), nullable=False)  # Amplifier Readout Noise
    naxis1 = Column(SmallInteger, nullable=False)  # Axis 1 size (pix)
    naxis2 = Column(SmallInteger, nullable=False)  # Axis 2 size (pix)
    ra_min = Column(Float(53), nullable=True)  # Image corner ra min (deg)
    ra_max = Column(Float(53), nullable=True)  # Image corner ra max (deg)
    dec_min = Column(Float(53), nullable=True)  # Image corner dec min (deg)
    dec_max = Column(Float(53), nullable=True)  # Image corner dec min (deg)
    zp_nightly = Column(Float(24), nullable=True)  # Zeropoint magnitude computed at the Nightly pipeline
    zp_nightly_err = Column(Float(24), nullable=True)  # Zeropoint magnitude error computed at the Nightly pipeline
    psf_fwhm = Column(Float(24), nullable=True)  # PSF FWHM measured on image. Available in reduced image only.
    cqa_1 = Column(Float(24), nullable=True)  # Quality Analysis check at observatory 1
    cqa_2 = Column(Float(24), nullable=True)  # Quality Analysis check at observatory 2
    cqa_3 = Column(Float(24), nullable=True)  # Quality Analysis check at observatory 3
    cqa_4 = Column(Float(24), nullable=True)  # Quality Analysis check at observatory 4
    cqa_5 = Column(Float(24), nullable=True)  # Quality Analysis check at observatory 5

    # Relationships
    mosaic = relationship('Mosaic', back_populates="images")
    detections = relationship('Detection', back_populates="image")
    detections_by_id = relationship('Detection', collection_class=attribute_mapped_collection('id'))
    forced_apertures = relationship('ForcedAperture', back_populates="image")

Index('ik_imagelocation', Image.ra_min, Image.ra_max, Image.dec_min, Image.dec_max)

    # Comment:
    # Contains the list of images associated to the mosaics (CCD and single amplifier images).


class Detection(Base):
    __tablename__ = 'detection'
    __table_args__ = (
        # Constraints
        PrimaryKeyConstraint('id'),
        ForeignKeyConstraint(['image_id'],          ['image.id'],      onupdate='CASCADE', ondelete='CASCADE'),
        )
    # Keys
    id = Column(BigInteger, nullable=False)  # Unique identifier
    image_id = Column(BigInteger, nullable=False)  # CCD image number
    insert_date = Column(DateTime, nullable=False, default=func.current_timestamp())  # Timestamp of insertion
    band = Column(String(8), nullable=False)  # Band name
    # Fields
    background = Column(Float(24), nullable=False)  # Background at centroid position [BACKGROUND] (count)
    class_star = Column(Float(24), nullable=False)  # Star-Galaxy classifier output [CLASS_STAR]
    spread_model = Column(Float(24), nullable=True)  # Spread parameter from model-fitting [SPREAD_MODEL]
    spreaderr_model = Column(Float(24), nullable=True)  # Spread parameter error from model-fitting [SPREADERR_MODEL]
    flux_auto = Column(Float(24), nullable=False)  # Flux within a Kron-like elliptical aperture [FLUX_AUTO] (count)
    flux_err_auto = Column(Float(24), nullable=False)  # RMS error for AUTO flux  [FLUXERR_AUTO] (count)
    flux_psf = Column(Float(24), nullable=False)  # Flux fitted using the PSF model (count)
    flux_err_psf = Column(Float(24), nullable=False)  # RMS error for PSF flux
    flux_model = Column(Float(24), nullable=False)  # Flux fitted using a galaxy model (count)
    flux_err_model = Column(Float(24), nullable=False)  # RMS error for the model flux
    flags = Column(SmallInteger, nullable=False)  # Extraction flags [FLAGS]
    elongation = Column(Float(24), nullable=False)  # A_IMAGE/B_IMAGE  [ELONGATION]
    dec = Column(Float(53), nullable=False)  # Windowed declination of barycenter (J2000) [DELTAWIN_J2000] (deg)
    ra = Column(Float(53), nullable=False)  # Windowed right ascension of barycenter (J2000) [ALPHAWIN_J2000] (deg)
    x = Column(Float(24), nullable=False)  # Windowed position estimate along x [XWIN_IMAGE] (pix)
    y = Column(Float(24), nullable=False)  # Windowed position estimate along y [YWIN_IMAGE] (pix)
    zp_offset = Column(Float(24), nullable=False)  # Offset zeropoint magnitude for unexpected corrections
    # Relationships
    image = relationship('Image', back_populates="detections")
    global_objects = relationship('Global_object',
                                  back_populates="detections",
                                  secondary=lambda: Global_object_detections.__table__)

Index('ik_detlocation', Detection.ra, Detection.dec)

    # Comment:
    # Contains the detections measured directly on the image after the nightly data reduction.

class ForcedAperture(Base):
    __tablename__ = 'forced_aperture'
    __table_args__ = (
        # Constraints
        PrimaryKeyConstraint('id'),
        ForeignKeyConstraint(['production_id'], ['production.id'], onupdate='CASCADE', ondelete='CASCADE'),
        ForeignKeyConstraint(['image_id'], ['image.id'], onupdate='CASCADE', ondelete='CASCADE'),
        )
    # Keys
    id = Column(BigInteger, nullable=False)  # Unique identifier
    production_id = Column(Integer, nullable=False)  # Production number
    pixel_id = Column(Integer, nullable=False)  # Healpix id (for area selection)
    ref_id = Column(BigInteger, nullable=True)  # Unique identifier for the reference object
    ref_cat = Column(String(16), nullable=True)  # Reference catalogue name
    image_id = Column(BigInteger, nullable=False)  # CCD image number
    insert_date = Column(DateTime, nullable=False, default=func.current_timestamp())  # Timestamp of insertion
    band = Column(String(8), nullable=False)  # Band name
    # Fields
    aperture_ra = Column(Float(53), nullable=False)  # Sky coordinates of aperture center
    aperture_dec = Column(Float(53), nullable=False)  # Sky coordinates of aperture center
    aperture_x = Column(Float(24), nullable=False)  # Image coordinates of aperture center
    aperture_y = Column(Float(24), nullable=False)  # Image coordinates of aperture center
    source_intensity = Column(Float(24), nullable=False)  # Integrated intensity of the source (in the aperture)
    source_uncertainty = Column(Float(24), nullable=False)  # Uncertainty associated with the source intensity
    magnitude = Column(Float(24), nullable=False)  # Magnitude representation of “SourceIntensity”.
    mag_uncertainty = Column(Float(24), nullable=False)  # Magnitude uncertainty, given by 1.0857 times
                                                         # “SourceUncertainty” divided by “SourceIntensity”.
    sky_median = Column(Float(24), nullable=False)  # The per-pixel median of samples in the sky annulus after the sky outliers have been rejected
    sky_sigma = Column(Float(24), nullable=False)  # The standard deviation of samples in the sky annulus after the sky outliers have been rejected
    radial_profile_fwhm = Column(Float(24), nullable=False)  # Full width at half maximum (FWHM) of the radial profile of the source (pixels).
    flag = Column(Boolean, nullable=False)  # Flag from pixel mask (False = Aperture OK, True = Aperture KO)
    # Relationships
    image = relationship('Image', back_populates="forced_apertures")
    production = relationship('Production', back_populates="forced_apertures")
    # Documentation

Index('ik_falocation', ForcedAperture.aperture_ra, ForcedAperture.aperture_dec)
Index('ik_fapixels', ForcedAperture.pixel_id)
Index('ik_faband', ForcedAperture.band)

    # Comment:
    # Contains the individual image measurements using force photometry in MEMBA.


# Forced Aperture Coadd Object Table
class ForcedApertureCoadd(Base):
    __tablename__ = 'forced_aperture_coadd'
    __table_args__ = (
        # Constraints
        PrimaryKeyConstraint('id'),
        ForeignKeyConstraint(['production_id'], ['production.id'], onupdate='CASCADE', ondelete='CASCADE'),
        UniqueConstraint('production_id', 'ra', 'dec'),
            )
    # Keys
    id = Column(                    BigInteger,   nullable=False) #Unique identifier"),
    production_id = Column(         Integer,      nullable=False) #Production number"),
    # Fields
    ref_id = Column(                    BigInteger,   nullable=False) #Unique identifier for reference catalogue"),
    ref_cat = Column(String(16), nullable=True)  # Reference catalogue name
    ra = Column(                    Float(53),    nullable=False) #Windowed right ascension of barycenter (J2000) [ALPHAWIN_J2000] (deg)"),
    dec = Column(                   Float(53),    nullable=False) #Windowed declination of barycenter (J2000) [DELTAWIN_J2000] (deg)"),
    # Magnitudes
    mag_u = Column(           Float(24),    nullable=True) #Magnitude fitted using a galaxy model"),
    mag_err_u = Column(       Float(24),    nullable=True) #RMS error for the model Magnitude"),
    mag_g = Column(           Float(24),    nullable=True) #Magnitude fitted using a galaxy model"),
    mag_err_g = Column(       Float(24),    nullable=True) #RMS error for the model Magnitude"),
    mag_r = Column(           Float(24),    nullable=True) #Magnitude fitted using a galaxy model"),
    mag_err_r = Column(       Float(24),    nullable=True) #RMS error for the model Magnitude"),
    mag_i = Column(           Float(24),    nullable=True) #Magnitude fitted using a galaxy model"),
    mag_err_i = Column(       Float(24),    nullable=True) #RMS error for the model Magnitude"),
    mag_z = Column(           Float(24),    nullable=True) #Magnitude fitted using a galaxy model"),
    mag_err_z = Column(       Float(24),    nullable=True) #RMS error for the model Magnitude"),
    mag_Y = Column(           Float(24),    nullable=True) #Magnitude fitted using a galaxy model"),
    mag_err_Y = Column(       Float(24),    nullable=True) #RMS error for the model Magnitude"),
    mag_NB455 = Column(         Float(24),    nullable=True) #Magnitude fitted using a galaxy model"),
    mag_err_NB455 = Column(     Float(24),    nullable=True) #RMS error for the model Magnitude"),
    mag_NB465 = Column(         Float(24),    nullable=True) #Magnitude fitted using a galaxy model"),
    mag_err_NB465 = Column(     Float(24),    nullable=True) #RMS error for the model Magnitude"),
    mag_NB475 = Column(         Float(24),    nullable=True) #Magnitude fitted using a galaxy model"),
    mag_err_NB475 = Column(     Float(24),    nullable=True) #RMS error for the model Magnitude"),
    mag_NB485 = Column(         Float(24),    nullable=True) #Magnitude fitted using a galaxy model"),
    mag_err_NB485 = Column(     Float(24),    nullable=True) #RMS error for the model Magnitude"),
    mag_NB495 = Column(         Float(24),    nullable=True) #Magnitude fitted using a galaxy model"),
    mag_err_NB495 = Column(     Float(24),    nullable=True) #RMS error for the model Magnitude"),
    mag_NB505 = Column(         Float(24),    nullable=True) #Magnitude fitted using a galaxy model"),
    mag_err_NB505 = Column(     Float(24),    nullable=True) #RMS error for the model Magnitude"),
    mag_NB515 = Column(         Float(24),    nullable=True) #Magnitude fitted using a galaxy model"),
    mag_err_NB515 = Column(     Float(24),    nullable=True) #RMS error for the model Magnitude"),
    mag_NB525 = Column(         Float(24),    nullable=True) #Magnitude fitted using a galaxy model"),
    mag_err_NB525 = Column(     Float(24),    nullable=True) #RMS error for the model Magnitude"),
    mag_NB535 = Column(         Float(24),    nullable=True) #Magnitude fitted using a galaxy model"),
    mag_err_NB535 = Column(     Float(24),    nullable=True) #RMS error for the model Magnitude"),
    mag_NB545 = Column(         Float(24),    nullable=True) #Magnitude fitted using a galaxy model"),
    mag_err_NB545 = Column(     Float(24),    nullable=True) #RMS error for the model Magnitude"),
    mag_NB555 = Column(         Float(24),    nullable=True) #Magnitude fitted using a galaxy model"),
    mag_err_NB555 = Column(     Float(24),    nullable=True) #RMS error for the model Magnitude"),
    mag_NB565 = Column(         Float(24),    nullable=True) #Magnitude fitted using a galaxy model"),
    mag_err_NB565 = Column(     Float(24),    nullable=True) #RMS error for the model Magnitude"),
    mag_NB575 = Column(         Float(24),    nullable=True) #Magnitude fitted using a galaxy model"),
    mag_err_NB575 = Column(     Float(24),    nullable=True) #RMS error for the model Magnitude"),
    mag_NB585 = Column(         Float(24),    nullable=True) #Magnitude fitted using a galaxy model"),
    mag_err_NB585 = Column(     Float(24),    nullable=True) #RMS error for the model Magnitude"),
    mag_NB595 = Column(         Float(24),    nullable=True) #Magnitude fitted using a galaxy model"),
    mag_err_NB595 = Column(     Float(24),    nullable=True) #RMS error for the model Magnitude"),
    mag_NB605 = Column(         Float(24),    nullable=True) #Magnitude fitted using a galaxy model"),
    mag_err_NB605 = Column(     Float(24),    nullable=True) #RMS error for the model Magnitude"),
    mag_NB615 = Column(         Float(24),    nullable=True) #Magnitude fitted using a galaxy model"),
    mag_err_NB615 = Column(     Float(24),    nullable=True) #RMS error for the model Magnitude"),
    mag_NB625 = Column(         Float(24),    nullable=True) #Magnitude fitted using a galaxy model"),
    mag_err_NB625 = Column(     Float(24),    nullable=True) #RMS error for the model Magnitude"),
    mag_NB635 = Column(         Float(24),    nullable=True) #Magnitude fitted using a galaxy model"),
    mag_err_NB635 = Column(     Float(24),    nullable=True) #RMS error for the model Magnitude"),
    mag_NB645 = Column(         Float(24),    nullable=True) #Magnitude fitted using a galaxy model"),
    mag_err_NB645 = Column(     Float(24),    nullable=True) #RMS error for the model Magnitude"),
    mag_NB655 = Column(         Float(24),    nullable=True) #Magnitude fitted using a galaxy model"),
    mag_err_NB655 = Column(     Float(24),    nullable=True) #RMS error for the model Magnitude"),
    mag_NB665 = Column(         Float(24),    nullable=True) #Magnitude fitted using a galaxy model"),
    mag_err_NB665 = Column(     Float(24),    nullable=True) #RMS error for the model Magnitude"),
    mag_NB675 = Column(         Float(24),    nullable=True) #Magnitude fitted using a galaxy model"),
    mag_err_NB675 = Column(     Float(24),    nullable=True) #RMS error for the model Magnitude"),
    mag_NB685 = Column(         Float(24),    nullable=True) #Magnitude fitted using a galaxy model"),
    mag_err_NB685 = Column(     Float(24),    nullable=True) #RMS error for the model Magnitude"),
    mag_NB695 = Column(         Float(24),    nullable=True) #Magnitude fitted using a galaxy model"),
    mag_err_NB695 = Column(     Float(24),    nullable=True) #RMS error for the model Magnitude"),
    mag_NB705 = Column(         Float(24),    nullable=True) #Magnitude fitted using a galaxy model"),
    mag_err_NB705 = Column(     Float(24),    nullable=True) #RMS error for the model Magnitude"),
    mag_NB715 = Column(         Float(24),    nullable=True) #Magnitude fitted using a galaxy model"),
    mag_err_NB715 = Column(     Float(24),    nullable=True) #RMS error for the model Magnitude"),
    mag_NB725 = Column(         Float(24),    nullable=True) #Magnitude fitted using a galaxy model"),
    mag_err_NB725 = Column(     Float(24),    nullable=True) #RMS error for the model Magnitude"),
    mag_NB735 = Column(         Float(24),    nullable=True) #Magnitude fitted using a galaxy model"),
    mag_err_NB735 = Column(     Float(24),    nullable=True) #RMS error for the model Magnitude"),
    mag_NB745 = Column(         Float(24),    nullable=True) #Magnitude fitted using a galaxy model"),
    mag_err_NB745 = Column(     Float(24),    nullable=True) #RMS error for the model Magnitude"),
    mag_NB755 = Column(         Float(24),    nullable=True) #Magnitude fitted using a galaxy model"),
    mag_err_NB755 = Column(     Float(24),    nullable=True) #RMS error for the model Magnitude"),
    mag_NB765 = Column(         Float(24),    nullable=True) #Magnitude fitted using a galaxy model"),
    mag_err_NB765 = Column(     Float(24),    nullable=True) #RMS error for the model Magnitude"),
    mag_NB775 = Column(         Float(24),    nullable=True) #Magnitude fitted using a galaxy model"),
    mag_err_NB775 = Column(     Float(24),    nullable=True) #RMS error for the model Magnitude"),
    mag_NB785 = Column(         Float(24),    nullable=True) #Magnitude fitted using a galaxy model"),
    mag_err_NB785 = Column(     Float(24),    nullable=True) #RMS error for the model Magnitude"),
    mag_NB795 = Column(         Float(24),    nullable=True) #Magnitude fitted using a galaxy model"),
    mag_err_NB795 = Column(     Float(24),    nullable=True) #RMS error for the model Magnitude"),
    mag_NB805 = Column(         Float(24),    nullable=True) #Magnitude fitted using a galaxy model"),
    mag_err_NB805 = Column(     Float(24),    nullable=True) #RMS error for the model Magnitude"),
    mag_NB815 = Column(         Float(24),    nullable=True) #Magnitude fitted using a galaxy model"),
    mag_err_NB815 = Column(     Float(24),    nullable=True) #RMS error for the model Magnitude"),
    mag_NB825 = Column(         Float(24),    nullable=True) #Magnitude fitted using a galaxy model"),
    mag_err_NB825 = Column(     Float(24),    nullable=True) #RMS error for the model Magnitude"),
    mag_NB835 = Column(         Float(24),    nullable=True) #Magnitude fitted using a galaxy model"),
    mag_err_NB835 = Column(     Float(24),    nullable=True) #RMS error for the model Magnitude"),
    mag_NB845 = Column(         Float(24),    nullable=True) #Magnitude fitted using a galaxy model"),
    mag_err_NB845 = Column(     Float(24),    nullable=True) #RMS error for the model Magnitude"),
    # Fluxes
    flux_u = Column(           Float(24),    nullable=True) #fluxnitude fitted using a galaxy model"),
    flux_err_u = Column(       Float(24),    nullable=True) #RMS error for the model fluxnitude"),
    flux_g = Column(           Float(24),    nullable=True) #fluxnitude fitted using a galaxy model"),
    flux_err_g = Column(       Float(24),    nullable=True) #RMS error for the model fluxnitude"),
    flux_r = Column(           Float(24),    nullable=True) #fluxnitude fitted using a galaxy model"),
    flux_err_r = Column(       Float(24),    nullable=True) #RMS error for the model fluxnitude"),
    flux_i = Column(           Float(24),    nullable=True) #fluxnitude fitted using a galaxy model"),
    flux_err_i = Column(       Float(24),    nullable=True) #RMS error for the model fluxnitude"),
    flux_z = Column(           Float(24),    nullable=True) #fluxnitude fitted using a galaxy model"),
    flux_err_z = Column(       Float(24),    nullable=True) #RMS error for the model fluxnitude"),
    flux_Y = Column(           Float(24),    nullable=True) #fluxnitude fitted using a galaxy model"),
    flux_err_Y = Column(       Float(24),    nullable=True) #RMS error for the model fluxnitude"),
    flux_NB455 = Column(         Float(24),    nullable=True) #fluxnitude fitted using a galaxy model"),
    flux_err_NB455 = Column(     Float(24),    nullable=True) #RMS error for the model fluxnitude"),
    flux_NB465 = Column(         Float(24),    nullable=True) #fluxnitude fitted using a galaxy model"),
    flux_err_NB465 = Column(     Float(24),    nullable=True) #RMS error for the model fluxnitude"),
    flux_NB475 = Column(         Float(24),    nullable=True) #fluxnitude fitted using a galaxy model"),
    flux_err_NB475 = Column(     Float(24),    nullable=True) #RMS error for the model fluxnitude"),
    flux_NB485 = Column(         Float(24),    nullable=True) #fluxnitude fitted using a galaxy model"),
    flux_err_NB485 = Column(     Float(24),    nullable=True) #RMS error for the model fluxnitude"),
    flux_NB495 = Column(         Float(24),    nullable=True) #fluxnitude fitted using a galaxy model"),
    flux_err_NB495 = Column(     Float(24),    nullable=True) #RMS error for the model fluxnitude"),
    flux_NB505 = Column(         Float(24),    nullable=True) #fluxnitude fitted using a galaxy model"),
    flux_err_NB505 = Column(     Float(24),    nullable=True) #RMS error for the model fluxnitude"),
    flux_NB515 = Column(         Float(24),    nullable=True) #fluxnitude fitted using a galaxy model"),
    flux_err_NB515 = Column(     Float(24),    nullable=True) #RMS error for the model fluxnitude"),
    flux_NB525 = Column(         Float(24),    nullable=True) #fluxnitude fitted using a galaxy model"),
    flux_err_NB525 = Column(     Float(24),    nullable=True) #RMS error for the model fluxnitude"),
    flux_NB535 = Column(         Float(24),    nullable=True) #fluxnitude fitted using a galaxy model"),
    flux_err_NB535 = Column(     Float(24),    nullable=True) #RMS error for the model fluxnitude"),
    flux_NB545 = Column(         Float(24),    nullable=True) #fluxnitude fitted using a galaxy model"),
    flux_err_NB545 = Column(     Float(24),    nullable=True) #RMS error for the model fluxnitude"),
    flux_NB555 = Column(         Float(24),    nullable=True) #fluxnitude fitted using a galaxy model"),
    flux_err_NB555 = Column(     Float(24),    nullable=True) #RMS error for the model fluxnitude"),
    flux_NB565 = Column(         Float(24),    nullable=True) #fluxnitude fitted using a galaxy model"),
    flux_err_NB565 = Column(     Float(24),    nullable=True) #RMS error for the model fluxnitude"),
    flux_NB575 = Column(         Float(24),    nullable=True) #fluxnitude fitted using a galaxy model"),
    flux_err_NB575 = Column(     Float(24),    nullable=True) #RMS error for the model fluxnitude"),
    flux_NB585 = Column(         Float(24),    nullable=True) #fluxnitude fitted using a galaxy model"),
    flux_err_NB585 = Column(     Float(24),    nullable=True) #RMS error for the model fluxnitude"),
    flux_NB595 = Column(         Float(24),    nullable=True) #fluxnitude fitted using a galaxy model"),
    flux_err_NB595 = Column(     Float(24),    nullable=True) #RMS error for the model fluxnitude"),
    flux_NB605 = Column(         Float(24),    nullable=True) #fluxnitude fitted using a galaxy model"),
    flux_err_NB605 = Column(     Float(24),    nullable=True) #RMS error for the model fluxnitude"),
    flux_NB615 = Column(         Float(24),    nullable=True) #fluxnitude fitted using a galaxy model"),
    flux_err_NB615 = Column(     Float(24),    nullable=True) #RMS error for the model fluxnitude"),
    flux_NB625 = Column(         Float(24),    nullable=True) #fluxnitude fitted using a galaxy model"),
    flux_err_NB625 = Column(     Float(24),    nullable=True) #RMS error for the model fluxnitude"),
    flux_NB635 = Column(         Float(24),    nullable=True) #fluxnitude fitted using a galaxy model"),
    flux_err_NB635 = Column(     Float(24),    nullable=True) #RMS error for the model fluxnitude"),
    flux_NB645 = Column(         Float(24),    nullable=True) #fluxnitude fitted using a galaxy model"),
    flux_err_NB645 = Column(     Float(24),    nullable=True) #RMS error for the model fluxnitude"),
    flux_NB655 = Column(         Float(24),    nullable=True) #fluxnitude fitted using a galaxy model"),
    flux_err_NB655 = Column(     Float(24),    nullable=True) #RMS error for the model fluxnitude"),
    flux_NB665 = Column(         Float(24),    nullable=True) #fluxnitude fitted using a galaxy model"),
    flux_err_NB665 = Column(     Float(24),    nullable=True) #RMS error for the model fluxnitude"),
    flux_NB675 = Column(         Float(24),    nullable=True) #fluxnitude fitted using a galaxy model"),
    flux_err_NB675 = Column(     Float(24),    nullable=True) #RMS error for the model fluxnitude"),
    flux_NB685 = Column(         Float(24),    nullable=True) #fluxnitude fitted using a galaxy model"),
    flux_err_NB685 = Column(     Float(24),    nullable=True) #RMS error for the model fluxnitude"),
    flux_NB695 = Column(         Float(24),    nullable=True) #fluxnitude fitted using a galaxy model"),
    flux_err_NB695 = Column(     Float(24),    nullable=True) #RMS error for the model fluxnitude"),
    flux_NB705 = Column(         Float(24),    nullable=True) #fluxnitude fitted using a galaxy model"),
    flux_err_NB705 = Column(     Float(24),    nullable=True) #RMS error for the model fluxnitude"),
    flux_NB715 = Column(         Float(24),    nullable=True) #fluxnitude fitted using a galaxy model"),
    flux_err_NB715 = Column(     Float(24),    nullable=True) #RMS error for the model fluxnitude"),
    flux_NB725 = Column(         Float(24),    nullable=True) #fluxnitude fitted using a galaxy model"),
    flux_err_NB725 = Column(     Float(24),    nullable=True) #RMS error for the model fluxnitude"),
    flux_NB735 = Column(         Float(24),    nullable=True) #fluxnitude fitted using a galaxy model"),
    flux_err_NB735 = Column(     Float(24),    nullable=True) #RMS error for the model fluxnitude"),
    flux_NB745 = Column(         Float(24),    nullable=True) #fluxnitude fitted using a galaxy model"),
    flux_err_NB745 = Column(     Float(24),    nullable=True) #RMS error for the model fluxnitude"),
    flux_NB755 = Column(         Float(24),    nullable=True) #fluxnitude fitted using a galaxy model"),
    flux_err_NB755 = Column(     Float(24),    nullable=True) #RMS error for the model fluxnitude"),
    flux_NB765 = Column(         Float(24),    nullable=True) #fluxnitude fitted using a galaxy model"),
    flux_err_NB765 = Column(     Float(24),    nullable=True) #RMS error for the model fluxnitude"),
    flux_NB775 = Column(         Float(24),    nullable=True) #fluxnitude fitted using a galaxy model"),
    flux_err_NB775 = Column(     Float(24),    nullable=True) #RMS error for the model fluxnitude"),
    flux_NB785 = Column(         Float(24),    nullable=True) #fluxnitude fitted using a galaxy model"),
    flux_err_NB785 = Column(     Float(24),    nullable=True) #RMS error for the model fluxnitude"),
    flux_NB795 = Column(         Float(24),    nullable=True) #fluxnitude fitted using a galaxy model"),
    flux_err_NB795 = Column(     Float(24),    nullable=True) #RMS error for the model fluxnitude"),
    flux_NB805 = Column(         Float(24),    nullable=True) #fluxnitude fitted using a galaxy model"),
    flux_err_NB805 = Column(     Float(24),    nullable=True) #RMS error for the model fluxnitude"),
    flux_NB815 = Column(         Float(24),    nullable=True) #fluxnitude fitted using a galaxy model"),
    flux_err_NB815 = Column(     Float(24),    nullable=True) #RMS error for the model fluxnitude"),
    flux_NB825 = Column(         Float(24),    nullable=True) #fluxnitude fitted using a galaxy model"),
    flux_err_NB825 = Column(     Float(24),    nullable=True) #RMS error for the model fluxnitude"),
    flux_NB835 = Column(         Float(24),    nullable=True) #fluxnitude fitted using a galaxy model"),
    flux_err_NB835 = Column(     Float(24),    nullable=True) #RMS error for the model fluxnitude"),
    flux_NB845 = Column(         Float(24),    nullable=True) #fluxnitude fitted using a galaxy model"),
    flux_err_NB845 = Column(     Float(24),    nullable=True) #RMS error for the model fluxnitude"),
    # Chi2
    chi2_u = Column(           Float(24),    nullable=True) #chi2nitude fitted using a galaxy model"),
    chi2_g = Column(           Float(24),    nullable=True) #chi2nitude fitted using a galaxy model"),
    chi2_r = Column(           Float(24),    nullable=True) #chi2nitude fitted using a galaxy model"),
    chi2_i = Column(           Float(24),    nullable=True) #chi2nitude fitted using a galaxy model"),
    chi2_z = Column(           Float(24),    nullable=True) #chi2nitude fitted using a galaxy model"),
    chi2_Y = Column(           Float(24),    nullable=True) #chi2nitude fitted using a galaxy model"),
    chi2_NB455 = Column(         Float(24),    nullable=True) #chi2nitude fitted using a galaxy model"),
    chi2_NB465 = Column(         Float(24),    nullable=True) #chi2nitude fitted using a galaxy model"),
    chi2_NB475 = Column(         Float(24),    nullable=True) #chi2nitude fitted using a galaxy model"),
    chi2_NB485 = Column(         Float(24),    nullable=True) #chi2nitude fitted using a galaxy model"),
    chi2_NB495 = Column(         Float(24),    nullable=True) #chi2nitude fitted using a galaxy model"),
    chi2_NB505 = Column(         Float(24),    nullable=True) #chi2nitude fitted using a galaxy model"),
    chi2_NB515 = Column(         Float(24),    nullable=True) #chi2nitude fitted using a galaxy model"),
    chi2_NB525 = Column(         Float(24),    nullable=True) #chi2nitude fitted using a galaxy model"),
    chi2_NB535 = Column(         Float(24),    nullable=True) #chi2nitude fitted using a galaxy model"),
    chi2_NB545 = Column(         Float(24),    nullable=True) #chi2nitude fitted using a galaxy model"),
    chi2_NB555 = Column(         Float(24),    nullable=True) #chi2nitude fitted using a galaxy model"),
    chi2_NB565 = Column(         Float(24),    nullable=True) #chi2nitude fitted using a galaxy model"),
    chi2_NB575 = Column(         Float(24),    nullable=True) #chi2nitude fitted using a galaxy model"),
    chi2_NB585 = Column(         Float(24),    nullable=True) #chi2nitude fitted using a galaxy model"),
    chi2_NB595 = Column(         Float(24),    nullable=True) #chi2nitude fitted using a galaxy model"),
    chi2_NB605 = Column(         Float(24),    nullable=True) #chi2nitude fitted using a galaxy model"),
    chi2_NB615 = Column(         Float(24),    nullable=True) #chi2nitude fitted using a galaxy model"),
    chi2_NB625 = Column(         Float(24),    nullable=True) #chi2nitude fitted using a galaxy model"),
    chi2_NB635 = Column(         Float(24),    nullable=True) #chi2nitude fitted using a galaxy model"),
    chi2_NB645 = Column(         Float(24),    nullable=True) #chi2nitude fitted using a galaxy model"),
    chi2_NB655 = Column(         Float(24),    nullable=True) #chi2nitude fitted using a galaxy model"),
    chi2_NB665 = Column(         Float(24),    nullable=True) #chi2nitude fitted using a galaxy model"),
    chi2_NB675 = Column(         Float(24),    nullable=True) #chi2nitude fitted using a galaxy model"),
    chi2_NB685 = Column(         Float(24),    nullable=True) #chi2nitude fitted using a galaxy model"),
    chi2_NB695 = Column(         Float(24),    nullable=True) #chi2nitude fitted using a galaxy model"),
    chi2_NB705 = Column(         Float(24),    nullable=True) #chi2nitude fitted using a galaxy model"),
    chi2_NB715 = Column(         Float(24),    nullable=True) #chi2nitude fitted using a galaxy model"),
    chi2_NB725 = Column(         Float(24),    nullable=True) #chi2nitude fitted using a galaxy model"),
    chi2_NB735 = Column(         Float(24),    nullable=True) #chi2nitude fitted using a galaxy model"),
    chi2_NB745 = Column(         Float(24),    nullable=True) #chi2nitude fitted using a galaxy model"),
    chi2_NB755 = Column(         Float(24),    nullable=True) #chi2nitude fitted using a galaxy model"),
    chi2_NB765 = Column(         Float(24),    nullable=True) #chi2nitude fitted using a galaxy model"),
    chi2_NB775 = Column(         Float(24),    nullable=True) #chi2nitude fitted using a galaxy model"),
    chi2_NB785 = Column(         Float(24),    nullable=True) #chi2nitude fitted using a galaxy model"),
    chi2_NB795 = Column(         Float(24),    nullable=True) #chi2nitude fitted using a galaxy model"),
    chi2_NB805 = Column(         Float(24),    nullable=True) #chi2nitude fitted using a galaxy model"),
    chi2_NB815 = Column(         Float(24),    nullable=True) #chi2nitude fitted using a galaxy model"),
    chi2_NB825 = Column(         Float(24),    nullable=True) #chi2nitude fitted using a galaxy model"),
    chi2_NB835 = Column(         Float(24),    nullable=True) #chi2nitude fitted using a galaxy model"),
    chi2_NB845 = Column(         Float(24),    nullable=True) #chi2nitude fitted using a galaxy model"),
    # N
    n_u = Column(           Float(24),    nullable=True) #fluxnitude fitted using a galaxy model"),
    n_g = Column(           Float(24),    nullable=True) #fluxnitude fitted using a galaxy model"),
    n_r = Column(           Float(24),    nullable=True) #fluxnitude fitted using a galaxy model"),
    n_i = Column(           Float(24),    nullable=True) #fluxnitude fitted using a galaxy model"),
    n_z = Column(           Float(24),    nullable=True) #fluxnitude fitted using a galaxy model"),
    n_Y = Column(           Float(24),    nullable=True) #fluxnitude fitted using a galaxy model"),
    n_NB455 = Column(         Float(24),    nullable=True) #fluxnitude fitted using a galaxy model"),
    n_NB465 = Column(         Float(24),    nullable=True) #fluxnitude fitted using a galaxy model"),
    n_NB475 = Column(         Float(24),    nullable=True) #fluxnitude fitted using a galaxy model"),
    n_NB485 = Column(         Float(24),    nullable=True) #fluxnitude fitted using a galaxy model"),
    n_NB495 = Column(         Float(24),    nullable=True) #fluxnitude fitted using a galaxy model"),
    n_NB505 = Column(         Float(24),    nullable=True) #fluxnitude fitted using a galaxy model"),
    n_NB515 = Column(         Float(24),    nullable=True) #fluxnitude fitted using a galaxy model"),
    n_NB525 = Column(         Float(24),    nullable=True) #fluxnitude fitted using a galaxy model"),
    n_NB535 = Column(         Float(24),    nullable=True) #fluxnitude fitted using a galaxy model"),
    n_NB545 = Column(         Float(24),    nullable=True) #fluxnitude fitted using a galaxy model"),
    n_NB555 = Column(         Float(24),    nullable=True) #fluxnitude fitted using a galaxy model"),
    n_NB565 = Column(         Float(24),    nullable=True) #fluxnitude fitted using a galaxy model"),
    n_NB575 = Column(         Float(24),    nullable=True) #fluxnitude fitted using a galaxy model"),
    n_NB585 = Column(         Float(24),    nullable=True) #fluxnitude fitted using a galaxy model"),
    n_NB595 = Column(         Float(24),    nullable=True) #fluxnitude fitted using a galaxy model"),
    n_NB605 = Column(         Float(24),    nullable=True) #fluxnitude fitted using a galaxy model"),
    n_NB615 = Column(         Float(24),    nullable=True) #fluxnitude fitted using a galaxy model"),
    n_NB625 = Column(         Float(24),    nullable=True) #fluxnitude fitted using a galaxy model"),
    n_NB635 = Column(         Float(24),    nullable=True) #fluxnitude fitted using a galaxy model"),
    n_NB645 = Column(         Float(24),    nullable=True) #fluxnitude fitted using a galaxy model"),
    n_NB655 = Column(         Float(24),    nullable=True) #fluxnitude fitted using a galaxy model"),
    n_NB665 = Column(         Float(24),    nullable=True) #fluxnitude fitted using a galaxy model"),
    n_NB675 = Column(         Float(24),    nullable=True) #fluxnitude fitted using a galaxy model"),
    n_NB685 = Column(         Float(24),    nullable=True) #fluxnitude fitted using a galaxy model"),
    n_NB695 = Column(         Float(24),    nullable=True) #fluxnitude fitted using a galaxy model"),
    n_NB705 = Column(         Float(24),    nullable=True) #fluxnitude fitted using a galaxy model"),
    n_NB715 = Column(         Float(24),    nullable=True) #fluxnitude fitted using a galaxy model"),
    n_NB725 = Column(         Float(24),    nullable=True) #fluxnitude fitted using a galaxy model"),
    n_NB735 = Column(         Float(24),    nullable=True) #fluxnitude fitted using a galaxy model"),
    n_NB745 = Column(         Float(24),    nullable=True) #fluxnitude fitted using a galaxy model"),
    n_NB755 = Column(         Float(24),    nullable=True) #fluxnitude fitted using a galaxy model"),
    n_NB765 = Column(         Float(24),    nullable=True) #fluxnitude fitted using a galaxy model"),
    n_NB775 = Column(         Float(24),    nullable=True) #fluxnitude fitted using a galaxy model"),
    n_NB785 = Column(         Float(24),    nullable=True) #fluxnitude fitted using a galaxy model"),
    n_NB795 = Column(         Float(24),    nullable=True) #fluxnitude fitted using a galaxy model"),
    n_NB805 = Column(         Float(24),    nullable=True) #fluxnitude fitted using a galaxy model"),
    n_NB815 = Column(         Float(24),    nullable=True) #fluxnitude fitted using a galaxy model"),
    n_NB825 = Column(         Float(24),    nullable=True) #fluxnitude fitted using a galaxy model"),
    n_NB835 = Column(         Float(24),    nullable=True) #fluxnitude fitted using a galaxy model"),
    n_NB845 = Column(         Float(24),    nullable=True) #fluxnitude fitted using a galaxy model"),
    flags = Column(SmallInteger, nullable=True) #Extraction flags [FLAGS]"),
    star_flag = Column(     Boolean,    nullable=True)

    #Reletionships
    production    = relationship('Production',       back_populates="forced_aperture_coadds")

Index('ik_forcedcoaddlocation', ForcedApertureCoadd.production_id, ForcedApertureCoadd.ra, ForcedApertureCoadd.dec)

    # Comment:
    # Contains the combined measurements using force photometry in MEMBA for all bands and passes for each reference source.


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

Index('ik_globallocation', Global_object.ra, Global_object.dec)

    # Comment:
    # Contains the list of unique objects from the multiple Detections table. (Currently unused)


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

    # Comment:
    # Many to Many intermediate table between Detection and Global_Object. (Currently unused)


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
    
    # Comment:
    # Contains unique coadd objects extracted from coadd image tiles. (Currently unused)

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

    # Comment:
    # Contains photometric redshift measurements from Coadd Forced Aperture measurements using the BCNz code.


Index('ik_photozredshift', Photo_z.production_id, Photo_z.z)
    

class SDSS_Star(Base):
    __tablename__ = 'sdss_star'
    __table_args__ = (
        # Constraints
        PrimaryKeyConstraint('objID'),
            )     
        
    # Keys
    objID = Column(BigInteger,   nullable=False) #Unique identifier. Same as Unique SDSS identifier composed from [skyVersion,rerun,run,camcol,field,obj]"),
    thingId = Column(Integer,   nullable=False) #Unique identifier. Same as Unique SDSS identifier composed from [skyVersion,rerun,run,camcol,field,obj]"),
    # Fields
    ra = Column(Float(53),    nullable=False) #Right Ascension of the object [RAJ2000] (deg)"),
    dec = Column(Float(53),    nullable=False) #Declination of the object [DEJ2000] (deg)"),
    raErr = Column(Float(24),    nullable=False) #Error in RA (arcsec)"),
    decErr = Column(Float(24),    nullable=False) #Error in DEC (arcsec)"),
    clean = Column(Float(24),    nullable=False) #Error in DEC (arcsec)"),
    psfMag_u = Column(Float(24),    nullable=False) #Model magnitude in u filter"),
    psfMag_g = Column(Float(24),    nullable=False) #Model magnitude in g filter"),
    psfMag_r = Column(Float(24),    nullable=False) #Model magnitude in r filter"),
    psfMag_i = Column(Float(24),    nullable=False) #Model magnitude in i filter"),
    psfMag_z = Column(Float(24),    nullable=False) #Model magnitude in z filter"),
    psfMagErr_u = Column(Float(24),    nullable=False) #Mean error on umag"),
    psfMagErr_g = Column(Float(24),    nullable=False) #Mean error on gmag"),
    psfMagErr_r = Column(Float(24),    nullable=False) #Mean error on rmag"),
    psfMagErr_i = Column(Float(24),    nullable=False) #Mean error on imag"),
    psfMagErr_z = Column(Float(24),    nullable=False) #Mean error on zmag"),
    # Constraints
    PrimaryKeyConstraint('objID'),

Index('ik_sdsslocation', SDSS_Star.ra, SDSS_Star.dec)

    # Comment:
    # External table from SDSS DR12 (Star view). Stars for simulation and calibration.

class SDSS_SpecPhoto(Base):
    __tablename__ = 'sdss_spec_photo'
    __table_args__ = (
        # Constraints
        PrimaryKeyConstraint('specObjID'),
            )

    # Keys
    objID = Column(BigInteger, nullable=False)
    specObjID = Column(BigInteger, nullable=False)
    mjd = Column(Integer, nullable=False)
    plate = Column(Integer, nullable=False)
    tile = Column(String(24), nullable=False)
    fiberID = Column(Integer, nullable=False)
    survey = Column(String(24), nullable=False)
    mode = Column(Integer, nullable=False)
    # Fields
    ra = Column(Float(53), nullable=False)  # Right Ascension of the object
    dec = Column(Float(53), nullable=False)  # Declination of the object
    z = Column(Float(24), nullable=False)
    zErr = Column(Float(24), nullable=False)
    zWarning = Column(Integer, nullable=False)
    _class = Column(String(32), nullable=False)  # GALAXY, STAR or QSO
    subClass = Column(String(32), nullable=False)
    fiberMag_u = Column(Float(24), nullable=False)
    fiberMag_g = Column(Float(24), nullable=False)
    fiberMag_r = Column(Float(24), nullable=False)
    fiberMag_i = Column(Float(24), nullable=False)
    fiberMag_z = Column(Float(24), nullable=False)
    fiberMagErr_u = Column(Float(24), nullable=False)
    fiberMagErr_g = Column(Float(24), nullable=False)
    fiberMagErr_r = Column(Float(24), nullable=False)
    fiberMagErr_i = Column(Float(24), nullable=False)
    fiberMagErr_z = Column(Float(24), nullable=False)
    modelMag_u = Column(Float(24), nullable=False)
    modelMag_g = Column(Float(24), nullable=False)
    modelMag_r = Column(Float(24), nullable=False)
    modelMag_i = Column(Float(24), nullable=False)
    modelMag_z = Column(Float(24), nullable=False)
    modelMagErr_u = Column(Float(24), nullable=False)
    modelMagErr_g = Column(Float(24), nullable=False)
    modelMagErr_r = Column(Float(24), nullable=False)
    modelMagErr_i = Column(Float(24), nullable=False)
    modelMagErr_z = Column(Float(24), nullable=False)
    extinction_u = Column(Float(24), nullable=False)
    extinction_g = Column(Float(24), nullable=False)
    extinction_r = Column(Float(24), nullable=False)
    extinction_i = Column(Float(24), nullable=False)
    extinction_z = Column(Float(24), nullable=False)

    # Constraints
    PrimaryKeyConstraint('specObjID'),

    # Comment:
    # External table from SDSS DR12 (Spec_Photo view). Sources with spectra for forced photometry and validation.

Index('ik_sdss_speclocation', SDSS_SpecPhoto.ra, SDSS_SpecPhoto.dec)

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

    # Comment:
    # External table from USNO. Bright stars for simulation and masking.

class DEEP2(Base):
    __tablename__ = 'deep2'
    __table_args__ = (
        # Constraints
        PrimaryKeyConstraint('objno'),
            )     
        
    # Keys
    objno = Column(        BigInteger,   nullable=False) #DEEP2 object number"),
    # Fields
    ra = Column(        Float(53),    nullable=False) #Right Ascension (in decimal degrees, J2000)"),
    dec = Column(       Float(53),    nullable=False) #Declination (in decimal degrees, J2000),
    magb = Column(      Float(24),    nullable=False) #CFHT B-band magnitude (AB) from Coil et al. 2004"),
    magr = Column(      Float(24),    nullable=False) #CFHT R-band magnitude (AB) from Coil et al. 2004"),
    magi = Column(      Float(24),    nullable=False) #CFHT I-band magnitude (AB) from Coil et al. 2004"),
    magberr = Column(      Float(24),    nullable=False) #B-band magnitude error"),
    magrerr = Column(      Float(24),    nullable=False) #R-band magnitude error"),
    magierr = Column(      Float(24),    nullable=False) #I-band magnitude error"),
    rg = Column(      Float(24),    nullable=False) #estimated R-band radius of object (sigma of Guassian fit in units of pixels --- 0.207”/pix)"),
    e2 = Column(      Float(24),    nullable=False) #ellipticity defined as E2 = (1 - b/a)"),
    pa = Column(      Float(24),    nullable=False) #object PA (degrees E of N)"),
    pgal = Column(      Float(24),    nullable=False) #the probability (0 - 1) that the sources is a galaxy for unresolved galaxies, 3 if resolved"),
    sfd_ebv = Column(      Float(24),    nullable=False) #E(B-V) from Schlegel, Finkbeiner, and Davis dust map"),
    m_b = Column(      Float(24),    nullable=False) #absolute B-band magnitude (AB, h = 1) from Willmer et al. (2006)"),
    ub = Column(      Float(24),    nullable=False) #rest-frame U-B color (AB) from Willmer et al. (2006)"),
    objname = Column(      String(8),    nullable=False) #the 8-digit DEEP2 object id (not always the same as OBJNO)"),
    mask = Column(      BigInteger,    nullable=False) #the DEEP2/DEEP3 slitmask number on which the object was observed"),
    slit = Column(      BigInteger,    nullable=False) #the slitlet number (on mask MASKNAME) in which the object was placed"),
    date = Column(      Date,    nullable=False) #Date on which the mask was observed (YYYY-MM-DD)"),
    mjd = Column(      Float(24),    nullable=False) #Modified Julian date of observation"),
    slitra = Column(      Float(24),    nullable=False) #RA of slit center"),
    slitdec = Column(      Float(24),    nullable=False) #Dec of slit center"),
    slitpa = Column(      Float(24),    nullable=False) #PA (degrees E of N) of slit"),
    slitlen = Column(      Float(24),    nullable=False) #length of slit (arcsec)"),
    z = Column(      Float(24),    nullable=False) #observed best-fitting redshift"),
    zbest = Column(      Float(24),    nullable=False) #best redshift (corrected for heliocentric motion)"),
    zerr = Column(      Float(24),    nullable=False) #redshift error (zerr < 0 indicates problematic z fit)"),
    zquality = Column(      Integer,    nullable=False) #redshift quality code, Q"),
    obj_type = Column(      String(6),    nullable=True) #type of best-fitting template (e.g., GALAXY or STAR)"),
    star_type = Column(      String(6),    nullable=True) #coarse classification for stellar templates"),
    rchi2 = Column(      Float(24),    nullable=False) #reduced chi-squared value for the redshift fit"),
    dof = Column(      BigInteger,    nullable=False) #degrees of freedom for redshift fit"),
    vdisp = Column(      Float(24),    nullable=False) #velocity dispersion in km/s"),
    vdisperr = Column(      Float(24),    nullable=False) #error in velocity dispersion"),
    comment = Column(      String(47),    nullable=True) #comment field"),
    # Constraints
    PrimaryKeyConstraint('objno'),

Index('ik_deeep2location', DEEP2.ra, DEEP2.dec)

    # Comment:
    # External table from DEEP2 Redshift Survey (DR4). Sources with spectra for forced photometry and validation.

class CFHTLenS(Base):
    __tablename__ = 'cfhtlens'
    __table_args__ = (
        # Constraints
        PrimaryKeyConstraint('paudm_id'),
            )

    # Keys
    paudm_id = Column(        BigInteger,   nullable=False) #PAUdm numeric id for reference table"),
    # Fields
    id = Column(        String(16),   nullable=False) #object id"),
    alpha_j2000 = Column(        Float(53),    nullable=False) #Right Ascension (in decimal degrees, J2000),
    delta_j2000 = Column(       Float(53),    nullable=False) #Declination (in decimal degrees, J2000),
    flag = Column(      Integer,    nullable=False) #whether the object has flags or not, 0: no flags, 1: some flags,
    class_star = Column(      Float(24),    nullable=False) #CLASS_STAR usually gives a cleaner star sample than star_flag, but can lead to serious incompleteness in a galaxy sample,
    e1 = Column(      Float(24),    nullable=False) #expectation values of galaxy ellipticity, from the marginalised galaxy ellipticity likelihood surface, to be used for shear measurement,
    e2 = Column(      Float(24),    nullable=False) #expectation values of galaxy ellipticity, from the marginalised galaxy ellipticity likelihood surface, to be used for shear measurement,
    scalelength = Column(      Float(24),    nullable=False) #lensfit galaxy model scalelength, marginalised over other parameters, as defined by Miller et al (2012),
    bulge_fraction = Column(      Float(24),    nullable=False) #ratio of the ﬂux in the bulge component to the total ﬂux (often written B/T)),
    snratio = Column(      Float(24),    nullable=False) #signal-to-noise ratio of the object, measured on a stack of the supplied exposures within a limiting isophote 2sigma above the noise,
    z_b = Column(      Float(24),    nullable=False) #most likely redshift from BPZ,
    z_b_min = Column(      Float(24),    nullable=False) #z_min 95% confidence interval from BPZ,
    z_b_max = Column(      Float(24),    nullable=False) #z_max 95% confidence interval from BPZ,
    t_b = Column(      Float(24),    nullable=False) #BPZ spectral type. 1=CWW-Ell, 2=CWW-Sbc, 3=CWW-Scd, 4=CWW-Im, 5=KIN-SB3, 6=KIN-SB2. Note that we use a recalibrated template set described in Capak et al. (2004) and that the templates are interpolated, hence fractional types occur,
    odds = Column(      Float(24),    nullable=False) #likelihood that z_b is correct from BPZ,
    z_ml = Column(      Float(24),    nullable=False) #maximum likelihood result (with "flat" unphysical prior) from BPZ,
    t_ml = Column(      Float(24),    nullable=False) #maximum likelihood result from BPZ,
    chi_squared_bpz = Column(      Float(24),    nullable=False) #"modified" chi square,
    star_flag = Column(      Integer,    nullable=False) #star_flag is optimized for galaxy studies, to keep an almost 100% complete galaxy sample with low (but not vanishing) stellar contamination,
    mag_u = Column(      Float(24),    nullable=False) #observed magnitude in the u-band,
    magerr_u = Column(      Float(24),    nullable=False) #error in the observed magnitude in the u-band,
    mag_g = Column(      Float(24),    nullable=False) #observed magnitude in the g-band,
    magerr_g = Column(      Float(24),    nullable=False) #error in the observed magnitude in the g-band,
    mag_r = Column(      Float(24),    nullable=False) #observed magnitude in the r-band,
    magerr_r = Column(      Float(24),    nullable=False) #error in the observed magnitude in the r-band,
    mag_i = Column(      Float(24),    nullable=False) #observed magnitude in the i-band,
    magerr_i = Column(      Float(24),    nullable=False) #error in the observed magnitude in the i-band,
    mag_y = Column(      Float(24),    nullable=False) #observed magnitude in the y-band,
    magerr_y = Column(      Float(24),    nullable=False) #error in the observed magnitude in the y-band,
    mag_z = Column(      Float(24),    nullable=False) #observed magnitude in the z-band,
    magerr_z = Column(      Float(24),    nullable=False) #error in the observed magnitude in the z-band,

Index('ik_cfhtlenslocation', CFHTLenS.alpha_j2000, CFHTLenS.delta_j2000)

    # Comment:
    # External table from CFHTLenS. Sources for forced photometry.


class COSMOS(Base):
    __tablename__ = 'cosmos'
    __table_args__ = (
        # Constraints
        PrimaryKeyConstraint('paudm_id'),
            )

    # Keys
    paudm_id = Column(BigInteger, nullable=False)  # PAUdm numeric id for reference table"),
    # Fields
    ra = Column(Float(53), nullable=False)  # (deg) Right Ascension in decimal degrees (J2000.0)
    dec = Column(Float(53), nullable=False)  # (deg) Declination in decimal degrees (J2000.0)
    zp_gal = Column(Float(24), nullable=False)  # photo-z using the galaxy templates (every source has a redshift)
    zl68_gal = Column(Float(24), nullable=False)  # lower limit, 68% confidence level
    zu68_gal = Column(Float(24), nullable=False)  # upper limit, 68% confidence level
    zl99_gal = Column(Float(24), nullable=False)  # lower li,mit, 99% confidence level
    zu99_gal = Column(Float(24), nullable=False)  # upper limit, 99% confidence level
    zp_sec = Column(Float(24), nullable=False)  # second photo-z solution
    dchi = Column(Float(24), nullable=False)  # chi_sec- chi_gal   reduced chi2  (-99 if less than 3 filters)
    Imag = Column(Float(24), nullable=False)  # Subaru i+   aperture magnitude (3" aperture)      7629.1 +- 1489.4/2
    eI = Column(Float(24), nullable=False)  # error on Imag  magnitudes
    I_auto = Column(Float(24), nullable=False)  # Auto magnitude  in the i band
    NbFilt = Column(Integer, nullable=False)  # Number of filters used in the photo-z fitting
    mod_gal = Column(Integer, nullable=False)  # best fit template  (1->8 Ell-S0, 9->15 Sa-Sc, 16->19 Sd-Sdm from Polletta et al., >=20 from BC03).
    type = Column(Integer, nullable=False)  # 0 if gal, 1 if star, 2 if Xray, 3 if IA>25.5, -99 if masked
    Umag = Column(Float(24), nullable=False)  # CFHT u*   aperture magnitude (3" aperture)        3911-0 +- 538/2
    Bmag = Column(Float(24), nullable=False)  # Subaru B   aperture magnitude (3" aperture)       4439.6 +-806.7/2
    Vmag = Column(Float(24), nullable=False)  # Subaru V   aperture magnitude (3" aperture)       5448.9 +- 934.8 /2
    Gmag = Column(Float(24), nullable=False)  # Subaru g+   aperture magnitude (3" aperture)      4728.3 +- 1162.9/2
    Rmag = Column(Float(24), nullable=False)  # Subaru r+   aperture magnitude (3" aperture)      6231.8 +- 1348.8/2
    Zmag = Column(Float(24), nullable=False)  # Subaru z+   aperture magnitude (3" aperture)      9021.6 +- 955.3/2
    ICmag = Column(Float(24), nullable=False)  # CFHT i'   aperture magnitude (3" aperture)       7628.9 +- 1460.0/2
    Jmag = Column(Float(24), nullable=False)  # UKIRT J   aperture magnitude (3" aperture)       12444.1 +- 1558.0/2
    Kmag = Column(Float(24), nullable=False)  # CFHT  K   aperture magnitude (3" aperture)       21480.2 +- 3250.0 /2
    MV = Column(Float(24), nullable=False)  # absolute magnitude in the V(Subaru) band
    ebv_gal = Column(Float(24), nullable=False)  # galactic extinction  E(B-V)
    ebv_int = Column(Float(24), nullable=False)  # additional extinction  E(B-V)
    zspec = Column(Float(24), nullable=False)  # spectroscopic redshift from ADP.2015-10-15T12_55_10.777_t2.txt
    conf = Column(Float(24), nullable=False)  # confidence class from ADP.2015-10-15T12_55_10.777_t2.txt
    F814W = Column(Float(24), nullable=False)  # HST  target magnitude from ADP.2015-10-15T12_55_10.777_t2.txt
    zfits = Column(String(60), nullable=True)  # filename of associated 1-d spectra ADP.2015-10-15T12_55_10.777_t2.txt

Index('ik_cosmoslocation', COSMOS.ra, COSMOS.dec)

    # Comment:
    # External table from zCOSMOS (DR3). Sources with accurate redshifts for forced photometry and validation.

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
Index('ik_truthlocation', Truth_Object.ra, Truth_Object.dec)

    # Comment:
    # (Simulation) Contains error-free simulated sources (stars and galaxies) for the pixel simulation pipeline.
    
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

Index('ik_targetlocation', Target.ra, Target.dec)

    # Comment:
    # (Simulation) Contains simulated exposure targets for the pixel simulation pipeline.

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

    # Comment:
    # (Simulation) Contains simulated background magnitude measurements for each target for the simulation pipeline.
    
    ###########    OPERATION Tables    ###########
    
    
    # Job Table
# class Job_pau(brownthrower.model.Job) :
#         
#         quality_controls = relationship('Quality_control',  back_populates="job")   
#     
#     # Documentation
#     #comment="Job",

# Job Comment:
# Tracks the list of Brownthrower computing jobs (Operation table).

# Dependency Comment:
# Tracks the dependency between Brownthrower jobs (Operation table).

# Tag Comment:
# Contains tags for Brownthrower jobs (Operation table).



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
    ref = Column(               String(32),   nullable=False) #QC Reference code"),
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
    
    # Comment:
    # Contains quality control entries measured during the data reduction process.

    
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
