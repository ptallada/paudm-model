#!/usr/bin/env python
# -*- coding: utf-8 -*-
import datetime
import optparse
import yaml

from sqlalchemy import ForeignKeyConstraint, Index, PrimaryKeyConstraint, UniqueConstraint
from sqlalchemy import create_engine
from sqlalchemy.types import BigInteger, Boolean, Date, Enum, Float, Integer, SmallInteger, String, Text, Time, DateTime
from sqlalchemy.ext.compiler import compiles
from sqlalchemy.orm import relationship, sessionmaker, scoped_session, synonym, contains_eager
from sqlalchemy.orm.collections import attribute_mapped_collection
from sqlalchemy.orm import backref
from sqlalchemy.sql.expression import func
from .base import Column
from sqlalchemy.ext.hybrid import hybrid_property
import brownthrower.model

Base = brownthrower.model.Base
metadata = Base.metadata
tables = {}
mappers = {}

_salt = 'djhferuyunniurnlp097jlknf8holanhgnhhrf'


class Project(Base):
    __tablename__ = 'project'
    __table_args__ = (
        # Primary key
        PrimaryKeyConstraint('id'),
        # Unique key
        UniqueConstraint('name'),
    )
    # Columns
    id = Column(Integer, nullable=False)
    name = Column(String(16), nullable=False)
    description = Column(String(128), nullable=True)
    contact_name = Column(String(64), nullable=False)
    contact_email = Column(String(64), nullable=False)
    created_at = Column(Date, nullable=False)

    # Relationships
    obs_sets = relationship('Obs_set', secondary='obs_set__project', back_populates='projects')


class Obs_set__Project(Base):
    __tablename__ = 'obs_set__project'
    __table_args__ = (
        # Constraints
        PrimaryKeyConstraint('project_id', 'obs_set_id'),
        # ForeignKeyConstraint(['red_img_id'] ,['red_image.id'],onupdate='CASCADE', ondelete='CASCADE'),
        ForeignKeyConstraint(
            ['project_id'],
            ['project.id'],
            onupdate='CASCADE', ondelete='CASCADE'
        ),
        ForeignKeyConstraint(
            ['obs_set_id'],
            ['obs_set.id'],
            onupdate='CASCADE', ondelete='CASCADE'
        ),
    )

    # Columns
    project_id = Column(BigInteger, nullable=False)
    obs_set_id = Column(BigInteger, nullable=False)


class Obs_set(Base):
    __tablename__ = 'obs_set'
    __table_args__ = (
        # Primary key
        PrimaryKeyConstraint('id'),
        # Unique key
        UniqueConstraint('obs_set', 'instrument'),
    )
    # Columns
    id = Column(Integer, nullable=False)  #
    instrument = Column(String(16), nullable=False)  # Name of the instrument used at the observation (i.e. PAUCam1.0)
    log = Column(String(128), nullable=True)  # Camera XML log file path
    rjd_start = Column(Float(53), nullable=False, index=True)  # First exposure observation time from current obs_set
    rjd_stop = Column(Float(53), nullable=False)  # Last exposure observation time from current obs_set
    night = Column(Date, nullable=False)  # Observation date
    notes = Column(Text, nullable=False)  # Observer's notes and comments
    operator = Column(String(128), nullable=False)  # Observer's name
    obs_set = Column(String(128), nullable=False)  # Observation Set identifier, from header

    # Relationships
    mosaics = relationship('Mosaic', back_populates="obs_set")
    projects = relationship('Project', secondary='obs_set__project', back_populates='obs_sets')

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
        ForeignKeyConstraint(['job_id'], ['job.id']),
    )
    # Columns
    id = Column(Integer, nullable=False)  #
    input_production_id = Column(Integer, nullable=True)  # input release id (from configuration)
    pipeline = Column(String(32), nullable=True)  # Pipeline Name [pixelsim, nightly, memba, analysis]
    release = Column(String(64), nullable=False)  # Major release name [TESTX, DRX]
    software_version = Column(String(32), nullable=True)  # Package version [DCX, vX.X]
    _comments = Column('comments', Text, nullable=True)
    job_id = Column(Integer, nullable=True)  # id of the job that generates the current production
    created = Column(DateTime, nullable=False, default=func.current_timestamp())  # Timestamp of insertion

    # Relationships
    mosaics = relationship('Mosaic', back_populates="production")
    forced_aperture_coadds = relationship('ForcedApertureCoadd', back_populates="production")
    forced_apertures = relationship('ForcedAperture', back_populates="production")
    photoz_bcnzs = relationship('Photoz_BCNz', back_populates="production")
    truth_objects = relationship('Truth_Object', back_populates="production")
    targets = relationship('Target', back_populates="production")
    phot_zps = relationship('PhotZP', back_populates="production")
    crosstalk_diffs = relationship('CrosstalkDiff', back_populates="production")
    crosstalk_ratios = relationship('CrosstalkRatio', back_populates="production")
    memba_ref_cats = relationship('MEMBARefCat', back_populates="production")
    parent = relationship('Production', back_populates='children',
                          primaryjoin='Production.input_production_id == Production.id',
                          remote_side='Production.id')
    children = relationship('Production', back_populates='parent',
                            primaryjoin='Production.input_production_id == Production.id',
                            cascade='all, delete-orphan', passive_deletes=True)
    job = relationship('Job',
                       primaryjoin='Production.job_id == Job.id',
                       foreign_keys='[ Production.job_id ]',
                       remote_side='[ Job.id ]',
                       uselist=False,
                       backref=backref("production", uselist=False))

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
            'timestamp': tnow,
            'comment': value
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
        ForeignKeyConstraint(['obs_set_id'], ['obs_set.id'], onupdate='CASCADE', ondelete='CASCADE'),
        # UniqueConstraint('archivepath', 'filename'),
        UniqueConstraint('production_id', 'obs_set_id', 'kind', 'exp_num'),

    )
    # Columns
    id = Column(BigInteger, nullable=False)  # Identifier
    production_id = Column(Integer, nullable=False)  # Production identifier
    obs_set_id = Column(Integer, nullable=False)  # obs_set number
    filename = Column(String(128), nullable=True)  # File name
    astrocat_filename = Column(String(128), nullable=True)  # File name for astro catalogue
    psfmodel_filename = Column(String(128), nullable=True)  # File name for psf model
    comment = Column(Text, nullable=True)
    archivepath = Column(String(128), nullable=True)  # Path in the archive
    kind = Column(
        Enum(
            'ARC', 'BIAS', 'DARK', 'FLASH', 'FLAT', 'FOCUS', 'GLANCE', 'PUPIL', 'SCIENCE',
            'SCRATCH', 'SKY', 'TARGET', 'MBIAS', 'MFLAT', 'RED_SCI', 'RED_WEIGHT', 'RED_MASK', 'STACKEDFOCUS', 'TEST',
            name='mosaic_kind'
        ),
        nullable=False)  # Mosaic image type
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
    extinction = Column(Float(24), nullable=True)  # Extinction value measured in reduction
    # for photometric calibration. Available in reduced image only.
    extinction_err = Column(Float(24), nullable=True)  # Extinction error value. Available in reduced image only.
    wind_spd = Column(Float(24), nullable=True)  # Wind speed (kph)
    wind_dir = Column(Float(24), nullable=True)  # Wind direction (deg)
    amb_temp = Column(Float(24), nullable=True)  # Ambient temperature (deg C)
    humidity = Column(Float(24), nullable=True)  # Ambient relative humidity (percent)
    pressure = Column(Float(24), nullable=True)  # Barometic pressure (mbar)
    astro_contrast = Column(Float(24), nullable=True)  # Astrometry contrast
    astro_chi2 = Column(Float(24), nullable=True)  # Astrometry chi2 fit to reference
    astro_nstars = Column(Float(24), nullable=True)  # Astrometry number of stars matched
    astro_nstars_highsn = Column(Float(24), nullable=True)  # Astrometry number of high SN stars matched
    astro_ref_sigma = Column(Float(24), nullable=True)  # Astrometry sigma offsets to reference
    astro_href_sigma = Column(Float(24), nullable=True)  # Astrometry sigma offsets to high SN reference
    astro_ref_cat = Column(String(32), nullable=False)  # Astrometry reference catalogue used
    merged_mosaics = Column(Integer, nullable=True)  # Number of mosaics merged to form the actual mosaic (for masters)
    mean_psf_fwhm = Column(Float(24), nullable=True)  # Mean PSF FWHM measured. Available in reduced image only,

    # detrend_status
    # 0: ok
    # 1: high read noise
    # 2: high saturated pixels
    # 3: high cosmics density
    # 4: multiple issues
    detrend_status = Column(SmallInteger, nullable=True)

    # astro_status
    # 0: ok
    # 1: all-sky reference
    # 2: low number of astro stars
    # 3: bad fit or low contrast
    # 4: no astrometry
    astro_status = Column(SmallInteger, nullable=True)

    # psf_model_status:
    # 0: ok
    # 1: low model stars
    # 2: detector PSF model failure
    # 3: focal plane PSF model failure
    psf_model_status = Column(SmallInteger, nullable=True)

    # photo_status:
    # 0: ok
    # 1: very high extinction
    # 2: detector photometry failure
    # 3: focal plane photometry failure
    # 4: no overlap with sdss
    photo_status = Column(SmallInteger, nullable=True)

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
    filename = Column(String(128), nullable=True)  # File name
    archivepath = Column(String(128), nullable=True)  # Path in the archive
    image_num = Column(SmallInteger, nullable=False)  # Extension number
    ccd_num = Column(SmallInteger, nullable=False)  # CCD Number"
    amp_num = Column(SmallInteger, nullable=False)  # Amplifier number (-1 for full CCD)
    filter = Column(String(8), nullable=True)  # Filter name"),
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
    zp_nightly_stars = Column(Float(24), nullable=True)  # Number of stars in calibration
    psf_fwhm = Column(Float(24), nullable=True)  # PSF FWHM measured on image. Available in reduced image only.
    bkg_mean = Column(Float(24), nullable=True)  # Mean background level (e-/s)
    bkg_std = Column(Float(24), nullable=True)  # STD background
    max_readnoise = Column(Float(24), nullable=True)  # Maximum readout noise from 4 amplifiers
    cosmic_ratio = Column(Float(24), nullable=True)  # Ratio of pixels flagged as cosmic rays
    saturate_ratio = Column(Float(24), nullable=True)  # Ratio of pixels flagged as saturated
    psf_stars = Column(Integer, nullable=True)  # Number of stars used to model the PSF
    psf_fit = Column(Float(24), nullable=True)  # Chi2 Fit of the PSF model
    n_extracted = Column(Integer, nullable=True)  # Number of sources extracted
    transparency = Column(Float(24), nullable=True)  # Transparency from photo_zp

    # Relationships
    mosaic = relationship('Mosaic', back_populates="images")
    detections = relationship('Detection', back_populates="image")
    detections_by_id = relationship('Detection', collection_class=attribute_mapped_collection('id'))
    forced_apertures = relationship('ForcedAperture', back_populates="image")
    star_photometries = relationship('StarPhotometry', back_populates="image")
    image_zps = relationship('ImageZP', back_populates="image")
    image_diff_origs = relationship('CrosstalkDiff', back_populates="image_orig", foreign_keys='[CrosstalkDiff.image_orig_id]')
    image_diff_dests = relationship('CrosstalkDiff', back_populates="image_dest", foreign_keys='[CrosstalkDiff.image_dest_id]')

Index('ik_imagelocation', Image.ra_min, Image.ra_max, Image.dec_min, Image.dec_max)


# Comment:
# Contains the list of images associated to the mosaics (CCD and single amplifier images).


class Detection(Base):
    __tablename__ = 'detection'
    __table_args__ = (
        # Constraints
        PrimaryKeyConstraint('id'),
        ForeignKeyConstraint(['image_id'], ['image.id'], onupdate='CASCADE', ondelete='CASCADE'),
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


Index('ik_detlocation', Detection.ra, Detection.dec)


# Comment:
# Contains the detections measured directly on the image after the nightly data reduction.


class StarZP(Base):
    __tablename__ = 'star_zp'
    __table_args__ = (
        # Constraints
        PrimaryKeyConstraint('id'),
        ForeignKeyConstraint(['star_photometry_id'], ['star_photometry.id'], onupdate='CASCADE', ondelete='CASCADE'),
    )
    # Keys
    id = Column(BigInteger, nullable=False)  # Unique identifier
    star_photometry_id = Column(BigInteger, nullable=False)

    # Fields
    zp = Column(Float(24), nullable=False)  # individual star zeropoint value
    zp_error = Column(Float(24), nullable=False)  # individual star zeropoint error
    chi2 = Column(Float(24), nullable=False)  # individual star fit chi2
    calib_method = Column(String(16), nullable=False)  # Calibration Method

    # Relationships
    star_photometry = relationship('StarPhotometry', back_populates="star_zps")

    # Comment:
    # Contains the individual zeropoint measurements for each star matched with the reference catalog
    # during the nightly photometry

Index('ik__star_zp__star_photometry_id', StarZP.star_photometry_id)


class StarPhotometry(Base):
    __tablename__ = 'star_photometry'
    __table_args__ = (
        # Constraints
        PrimaryKeyConstraint('id'),
        ForeignKeyConstraint(['image_id'], ['image.id'], onupdate='CASCADE', ondelete='CASCADE'),
        ForeignKeyConstraint(['phot_method_id'], ['phot_method.id'], onupdate='CASCADE', ondelete='CASCADE'),
        # UniqueConstraint('image_id', 'ref_id', 'phot_method_id'),
    )
    # Keys
    id = Column(BigInteger, nullable=False)  # Unique identifier
    image_id = Column(BigInteger, nullable=False)  # CCD image number

    # Non-relationable Keys
    ref_cat = Column(String(16), nullable=False)  # Reference catalogue name
    ref_id = Column(BigInteger, nullable=True)  # Reference catalogue star id

    # Fields
    x_image = Column(Float(24), nullable=False)  # x position in the image
    y_image = Column(Float(24), nullable=False)  # y position in the image
    flux = Column(Float(24), nullable=False)  # measured flux
    flux_err = Column(Float(24), nullable=False)  # measured flux error
    bg = Column(Float(24), nullable=False)  # measured background flux
    bg_err = Column(Float(24), nullable=False)  # measured background flux error
    flags = Column(SmallInteger, nullable=True)  # Extraction flags
    phot_method_id = Column(Integer, nullable=False)  # Photometry method id

    # Relationships
    image = relationship('Image', back_populates="star_photometries")
    star_zps = relationship('StarZP', back_populates="star_photometry")
    phot_method = relationship('PhotMethod', back_populates="star_photometries")


    # Comment:
    # Contains the individual photometry measurements for each star matched with the reference catalogue
    # during the nightly photometry

Index('ik__star_photometry__image_id', StarPhotometry.image_id)
Index('ik__star_photometry__phot_method_id', StarPhotometry.phot_method_id)


class PhotMethod(Base):
    __tablename__ = 'phot_method'
    __table_args__ = (
        # Constraints
        PrimaryKeyConstraint('id'),
    )
    # Keys
    id = Column(Integer, nullable=False)  # Unique identifier

    # Extraction method
    extraction_code = Column(String(16), nullable=False)  # Code used for extraction (i.e. sextractor/photutils)
    extraction_method = Column(String(16), nullable=False)  # Extraction method (i.e. APER/AUTO)
    extraction_parameter = Column(Float(24), nullable=True)  # Optional parameter for extraction

    # Background method
    background_method = Column(String(16), nullable=False)  # Background method (global, local, annulus)
    background_parameter = Column(Float(24), nullable=True)  # Optional parameter for extraction

    # Scatterlight method
    scatterlight_method = Column(String(16), nullable=False)  # Scatterlight correction method
    scatterlight_parameter = Column(Float(24), nullable=True)  # Optional parameter for Scatterlight correction

    _comments = Column('comments', Text, nullable=True)

    star_photometries = relationship('StarPhotometry', back_populates="phot_method")
    image_zps = relationship('ImageZP', back_populates="phot_method")

    # Comment:
    # Contains the information of the photometry method


class ImageZP(Base):
    __tablename__ = 'image_zp'
    __table_args__ = (
        # Constraints
        PrimaryKeyConstraint('id'),
        ForeignKeyConstraint(['image_id'], ['image.id'], onupdate='CASCADE', ondelete='CASCADE'),
        ForeignKeyConstraint(['phot_method_id'], ['phot_method.id'], onupdate='CASCADE', ondelete='CASCADE'),
        # UniqueConstraint('image_id', 'ref_id', 'phot_method_id'),
    )
    # Keys
    id = Column(BigInteger, nullable=False)  # Unique identifier
    image_id = Column(BigInteger, nullable=False)  # CCD image number

    # Fields
    zp = Column(Float(24), nullable=False)  # Image zeropoint value
    zp_error = Column(Float(24), nullable=False)  # Image zeropoint error
    n_stars = Column(Integer, nullable=True)  # Number of stars used to build image zp
    phot_method_id = Column(Integer, nullable=False)  # Photometry method
    calib_method = Column(String(16), nullable=False)  # Calibration method
    transparency = Column(Float(24), nullable=True)  # Image zeropoint value

    # Relationships
    image = relationship('Image', back_populates="image_zps")
    phot_method = relationship('PhotMethod', back_populates="image_zps")

    # Comment:
    # Contains the image zeropoint measurements for each photometry-calibration method

Index('ik__image_zp__image_id', ImageZP.image_id)
Index('ik__image_zp__phot_method_id', ImageZP.phot_method_id)


class PhotZP(Base):
    __tablename__ = 'phot_zp'
    __table_args__ = (
        # Constraints
        PrimaryKeyConstraint('id'),
        ForeignKeyConstraint(['production_id'], ['production.id'], onupdate='CASCADE', ondelete='CASCADE'),
        UniqueConstraint('production_id', 'band'),

    )
    # Keys
    id = Column(BigInteger, nullable=False)  # Unique identifier
    production_id = Column(Integer, nullable=False)  # Production identifier

    # Fields
    zp = Column(Float(24), nullable=False)  # x position in the image
    band = Column(String(8), nullable=False)  # Band name
    date = Column(DateTime, nullable=False, default=func.current_timestamp())  # Timestamp of insertion

    # Relationships
    production = relationship('Production', back_populates="phot_zps")

    # Comment:
    # Contains the photometric zeropoints for a given production


class CrosstalkDiff(Base):
    __tablename__ = 'crosstalk_diff'
    __table_args__ = (
        # Constraints
        PrimaryKeyConstraint('image_orig_id', 'image_dest_id', 'production_id'),
        ForeignKeyConstraint(['image_orig_id'], ['image.id'], onupdate='CASCADE', ondelete='CASCADE'),
        ForeignKeyConstraint(['image_dest_id'], ['image.id'], onupdate='CASCADE', ondelete='CASCADE'),
        ForeignKeyConstraint(['production_id'], ['production.id'], onupdate='CASCADE', ondelete='CASCADE'),
    )
    # Keys
    image_orig_id = Column(BigInteger, nullable=False)  # Origin image identifier
    image_dest_id = Column(BigInteger, nullable=False)  # Destination image indetifier
    production_id = Column(Integer, nullable=False)  # Production id

    # Fields
    background_orig_all = Column(Float(53), nullable=False)  # Background level at origin for all pixels
    background_dest_all = Column(Float(53), nullable=False)  # Background level at destination for all pixels
    background_dest_sat = Column(Float(53), nullable=False)  # Background level at destination for saturated pixels
    npix_orig_sat = Column(Integer, nullable=False)  # Number of pixels at origin

    # Relationships
    image_orig = relationship('Image', back_populates="image_diff_origs", foreign_keys=[image_orig_id])
    image_dest = relationship('Image', back_populates="image_diff_dests", foreign_keys=[image_dest_id])
    production = relationship('Production', back_populates="crosstalk_diffs")

Index('ik_crosstalkdiff', CrosstalkDiff.production_id, CrosstalkDiff.image_orig_id, CrosstalkDiff.image_dest_id)


    # Comment:
    # Crosstalk difference values


class CrosstalkRatio(Base):
    __tablename__ = 'crosstalk_ratio'
    __table_args__ = (
        # Constraints
        PrimaryKeyConstraint('ccd_num_orig', 'amp_num_orig', 'ccd_num_dest', 'amp_num_dest', 'production_id'),
        ForeignKeyConstraint(['production_id'], ['production.id'], onupdate='CASCADE', ondelete='CASCADE'),
    )
    # Keys
    ccd_num_orig = Column(SmallInteger, nullable=False)  # Origin CCD number
    amp_num_orig = Column(SmallInteger, nullable=False)  # Origin amplifier number
    ccd_num_dest = Column(SmallInteger, nullable=False)  # Destination CCD number
    amp_num_dest = Column(SmallInteger, nullable=False)  # Destination amplifier number
    production_id = Column(Integer, nullable=False)  # Production id

    # Fields
    ratio = Column(Float(53), nullable=False)  # crosstalk ratio between amplifiers

    # Relationships
    production = relationship('Production', back_populates="crosstalk_ratios")


Index('ik_crosstalkratio', CrosstalkRatio.production_id, CrosstalkRatio.ccd_num_orig, CrosstalkRatio.amp_num_orig,
                                                   CrosstalkRatio.ccd_num_dest, CrosstalkRatio.amp_num_dest)

    # Comment:
    # Crosstalk ratio values


# MEMBA production to Reference Catalog Table
class MEMBARefCat(Base):
    __tablename__ = 'memba_ref_cat'
    __table_args__ = (
        # Constraints
        PrimaryKeyConstraint('production_id', 'ref_cat'),
        ForeignKeyConstraint(['production_id'], ['production.id'], onupdate='CASCADE', ondelete='CASCADE'),
    )
    # Keys
    production_id = Column(Integer, nullable=False)  # Production id
    ref_cat = Column(String(16),  nullable=True)  # Reference catalogue name
    production = relationship('Production', back_populates="memba_ref_cats")


# Forced Aperture Object Table
class ForcedAperture(Base):
    __tablename__ = 'forced_aperture'
    __table_args__ = (
        # Constraints
        PrimaryKeyConstraint('production_id', 'ref_id', 'image_id'),
        ForeignKeyConstraint(['production_id'], ['production.id'], onupdate='CASCADE', ondelete='CASCADE'),
        ForeignKeyConstraint(['image_id'], ['image.id'], onupdate='CASCADE', ondelete='CASCADE'),
    )
    # Keys
    production_id = Column(Integer, nullable=False, comment='Production id')
    image_id = Column(BigInteger, nullable=False, comment='CCD image id')
    ref_id = Column(BigInteger, nullable=True, comment='Unique identifier of source in reference catalogue')
    # Fields
    pixel_id = Column(Integer, nullable=False, comment='Healpix id')
    aperture_x = Column(Float(24), nullable=False, comment='Image X coordinate of aperture center (pixels)')
    aperture_y = Column(Float(24), nullable=False, comment='Image Y coordinate of aperture center (pixels)')
    aperture_a = Column(Float(24), nullable=False, comment='aperture size radius for major axis (arcsec)')
    aperture_b = Column(Float(24), nullable=False, comment='aperture size radius for minor axis (arcsec)')
    aperture_theta = Column(Float(24), nullable=False, comment='aperture position angle from N to E(deg)')
    flux = Column(Float(24), nullable=False, comment='Non calibrated flux in the aperture (e/s)')
    flux_error = Column(Float(24), nullable=False, comment='Uncertainty associated to the non calibrated flux in the aperture')
    annulus_a_in = Column(Float(24), nullable=False, comment='annulus inner radius for major axis (arcsec)')
    annulus_a_out = Column(Float(24), nullable=False, comment='annulus outer radius for major axis (arcsec)')
    annulus_b_in = Column(Float(24), nullable=False, comment='annulus inner radius for minor axis (arcsec)')
    annulus_b_out = Column(Float(24), nullable=False, comment='annulus outer radius for minor axis (arcsec)')
    annulus_median = Column(Float(24), nullable=False, comment='Sky median in annulus after sigma clipping')
    annulus_sigma = Column(Float(24), nullable=False, comment='Sky standard deviation in annulus after sigma clipping')
    annulus_samples = Column(Integer, nullable=False, comment='The number of samples in the annulus used for sky statistics')
    annulus_ellipticity = Column(Float(24), nullable=True, comment='Local background ellipticity')
    image_ellipticity = Column(Float(24), nullable=True, comment='Global background ellipticity')
    flag = Column(Integer, nullable=True, comment='Flag value from mask and MEMBA analysis')

    # Relationships
    image = relationship('Image', back_populates="forced_apertures")
    production = relationship('Production', back_populates="forced_apertures")


# Indexes
Index('ik_faepixels', ForcedAperture.production_id, ForcedAperture.pixel_id)
Index('ik__fa__image_id', ForcedAperture.image_id)

# Comment
# Contains the individual image measurements using force photometry in MEMBA.


# Forced Aperture Coadd Object Table
class ForcedApertureCoadd(Base):
    __tablename__ = 'forced_aperture_coadd'
    __table_args__ = (
        # Constraints
        PrimaryKeyConstraint('production_id', 'ref_id', 'band'),
        ForeignKeyConstraint(['production_id'], ['production.id'], onupdate='CASCADE', ondelete='CASCADE'),
    )
    # Keys
    production_id = Column(Integer, nullable=False, comment='Production id')
    ref_id = Column(BigInteger, nullable=False, comment='Unique identifier of source in reference catalogue')
    band = Column(Enum('NB455', 'NB465', 'NB475', 'NB485', 'NB495',
                       'NB505', 'NB515', 'NB525', 'NB535', 'NB545', 'NB555', 'NB565', 'NB575', 'NB585', 'NB595',
                       'NB605', 'NB615', 'NB625', 'NB635', 'NB645', 'NB655', 'NB665', 'NB675', 'NB685', 'NB695',
                       'NB705', 'NB715', 'NB725', 'NB735', 'NB745', 'NB755', 'NB765', 'NB775', 'NB785', 'NB795',
                       'NB805', 'NB815', 'NB825', 'NB835', 'NB845', 'u', 'g', 'r', 'i', 'z', 'Y',
                       name='pau_bands'), nullable=False, comment='Band')
    # Fields
    flux = Column(Float(24), nullable=True, comment='Calibrated flux (e/s)')
    flux_error = Column(Float(24), nullable=True, comment='Calibrated flux error')
    chi2 = Column(Float(24), nullable=True, comment='Chi Square of fit from multiple observations')
    n_coadd = Column(SmallInteger, nullable=True, comment='Number of coadded observations')
    run = Column(SmallInteger, nullable=True, comment='run number for coadds in memba production')

    # Relationships
    production = relationship('Production', back_populates="forced_aperture_coadds")

# Indexes
# Already useful the primary keys


"""
# Forced Aperture Report table
class ForcedApertureReport(Base):
    __tablename__ = 'forced_aperture_report'
    __table_args__ = (
        # Constraints
        PrimaryKeyConstraint('id'),
        ForeignKeyConstraint(['fac_id'], ['forced_aperture_coadd.id'], onupdate='CASCADE', ondelete='CASCADE'),
    )

    # Keys
    id = Column(Integer, nullable=False)  # Unique identifier"),
    fac_id = Column(BigInteger, nullable=False)  # Forced aperture coadd id"),
    # Fields
    fa_id = Column(BigInteger, nullable=True)  # Forced aperture id (optional)
    band = Column(String(16), nullable=True)  # Reference catalogue name
    # REPORT status:
    # FAC global issues
    #  0: everything seems ok
    #  1: major global problem

    # FAC band issues
    #  100: flux underestimated
    #  101: flux overestimated
    #  102: noise underestimated
    #  103: noise overestimated
    #  104: discordant measurements

    # FA issues
    #  200: flag wrong
    #  201: flux underestimated
    #  202: flux overestimated
    #  203: noise underestimated
    #  204: noise overestimated
    #  205: missing bright source
    #  206: noisy image
    #  207: astrometry issue
    #  208: scatterlight
    #  209: edge effect
    #  210: bad pixel
    #  211: blended
    #  212: cosmic ray
    #  213: star halo/spikes
    #  214: zeropoint issue
    #  215: missing r50
    #  216: crosstalk
    #  217: other

    report_status = Column(Integer, nullable=False)  # Report id
    user = Column(Text, nullable=False)  # username
    insert_date = Column(DateTime, nullable=False, default=func.current_timestamp())  # Timestamp of insertion

    # Relationships
    fac = relationship('ForcedApertureCoadd', back_populates="forced_aperture_reports")


Index('ik_forcedreport', ForcedApertureReport.fac_id, ForcedApertureReport.band)
"""


class Photoz_BCNz(Base):
    __tablename__ = 'photoz_bcnz'
    __table_args__ = (
        # Constraints
        PrimaryKeyConstraint('production_id', 'ref_id'),
        ForeignKeyConstraint(['production_id'], ['production.id'], onupdate='CASCADE', ondelete='CASCADE'),
    )
    # Keys
    production_id = Column(Integer, nullable=False, comment='Production id')
    ref_id = Column(BigInteger, nullable=False, comment='Reference id')
    # Fields
    zb = Column(Float(24), nullable=False, comment='Bayesian photometric redshift (peak of p(z))')
    odds = Column(Float(24), nullable=False, comment='ODDS quality parameter')
    pz_width = Column(Float(24), nullable=False, comment='pz_width quality parameter')
    zb_mean = Column(Float(24), nullable=False, comment='Bayesian photometric redshift (mean of p(z))')
    chi2 = Column(Float(24), nullable=False, comment='Minimum chi2')
    n_band = Column(Float(24), nullable=False, comment='Number of bands used for chi2 fit')
    ebv = Column(Float(24), nullable=False, comment='E(B-V) extinction value')

    # Relationships
    production = relationship('Production', back_populates="photoz_bcnzs")

    # Comment:
    # Contains photometric redshift measurements using the BCNz v2 code.

Index('ik_photozredshift', Photoz_BCNz.production_id, Photoz_BCNz.zb)


class SDSS_Star(Base):
    __tablename__ = 'sdss_star'
    __table_args__ = (
        # Constraints
        PrimaryKeyConstraint('objID'),
    )

    # Keys
    objID = Column(BigInteger,
                   nullable=False)  # Unique identifier. Same as Unique SDSS identifier composed from [skyVersion,rerun,run,camcol,field,obj]"),
    thingId = Column(Integer,
                     nullable=False)  # Unique identifier. Same as Unique SDSS identifier composed from [skyVersion,rerun,run,camcol,field,obj]"),
    # Fields
    ra = Column(Float(53), nullable=False)  # Right Ascension of the object [RAJ2000] (deg)"),
    dec = Column(Float(53), nullable=False)  # Declination of the object [DEJ2000] (deg)"),
    raErr = Column(Float(24), nullable=False)  # Error in RA (arcsec)"),
    decErr = Column(Float(24), nullable=False)  # Error in DEC (arcsec)"),
    clean = Column(Float(24), nullable=False)  # Error in DEC (arcsec)"),
    psfMag_u = Column(Float(24), nullable=False)  # Model magnitude in u filter"),
    psfMag_g = Column(Float(24), nullable=False)  # Model magnitude in g filter"),
    psfMag_r = Column(Float(24), nullable=False)  # Model magnitude in r filter"),
    psfMag_i = Column(Float(24), nullable=False)  # Model magnitude in i filter"),
    psfMag_z = Column(Float(24), nullable=False)  # Model magnitude in z filter"),
    psfMagErr_u = Column(Float(24), nullable=False)  # Mean error on umag"),
    psfMagErr_g = Column(Float(24), nullable=False)  # Mean error on gmag"),
    psfMagErr_r = Column(Float(24), nullable=False)  # Mean error on rmag"),
    psfMagErr_i = Column(Float(24), nullable=False)  # Mean error on imag"),
    psfMagErr_z = Column(Float(24), nullable=False)  # Mean error on zmag"),
    # Constraints
    PrimaryKeyConstraint('objID'),


Index('ik_sdsslocation', SDSS_Star.ra, SDSS_Star.dec)


# Comment:
# External table from SDSS DR12 (Star view). Stars for simulation and calibration.


class gaia(Base):
    __tablename__ = 'gaia'
    __table_args__ = (
        # Constraints
        PrimaryKeyConstraint('source_id'),
            )     
        
    # Keys
    source_id = Column(BigInteger,   nullable=False) #Unique identifier.
    # Fields
    ra = Column(Float(53),    nullable=False) #Right Ascension of the object (deg)"),
    dec = Column(Float(53),    nullable=False) #Declination of the object (deg)"),
    ra_err = Column(Float(24),    nullable=False) #Error in RA (arcsec)"),
    dec_err = Column(Float(24),    nullable=False) #Error in DEC (arcsec)"),
    phot_g_mean_mag = Column(Float(24),    nullable=False) #Magnitude in g filter"),
    phot_g_mean_flux = Column(Float(24),    nullable=False) #Flux in g filter"),
    phot_g_mean_flux_error = Column(Float(24),    nullable=False) #Flux error in g filter"),
    ref_epoch = Column(Float(24),    nullable=False) #Ref epoch gaia"),
    # Constraints
    PrimaryKeyConstraint('source_id'),

Index('ik_gaialocation', gaia.ra, gaia.dec)

    # Comment:
    # External table from gaia. Stars for simulation and calibration.


class gaia_dr2(Base):
    __tablename__ = 'gaia_dr2'
    __table_args__ = (
        # Constraints
        PrimaryKeyConstraint('source_id'),
    )

    # Keys
    source_id = Column(BigInteger, nullable=False)  # Unique identifier.
    # Fields
    duplicated_source = Column(Boolean, nullable=False)
    ra = Column(Float(53), nullable=False)  # Right Ascension of the object (deg)
    dec = Column(Float(53), nullable=False)  # Declination of the object (deg)
    ra_err = Column(Float(24), nullable=False)  # Error in RA (arcsec)
    dec_err = Column(Float(24), nullable=False)  # Error in DEC (arcsec)
    pmra = Column(Float(24), nullable=False)  # Proper motion in RA
    pmdec = Column(Float(24), nullable=False)  # Proper motion in Dec
    pmra_error = Column(Float(24), nullable=False)  # Proper motion error in RA
    pmdec_error = Column(Float(24), nullable=False)  # Proper motion error in Dec
    phot_g_mean_mag = Column(Float(24), nullable=False)  # Magnitude in g filter
    phot_g_mean_flux = Column(Float(24), nullable=False)  # Flux in g filter
    phot_g_mean_flux_error = Column(Float(24), nullable=False)  # Flux error in g filter
    phot_bp_mean_mag = Column(Float(24), nullable=False)  # Magnitude in bp filter
    phot_bp_mean_flux = Column(Float(24), nullable=False)  # Flux in bp filter
    phot_bp_mean_flux_error = Column(Float(24), nullable=False)  # Flux error in bp filter
    phot_rp_mean_mag = Column(Float(24), nullable=False)  # Magnitude in rp filter
    phot_rp_mean_flux = Column(Float(24), nullable=False)  # Flux in rp filter
    phot_rp_mean_flux_error = Column(Float(24), nullable=False)  # Flux error in rp filter
    # Constraints
    PrimaryKeyConstraint('source_id'),


Index('ik_gaia_dr2location', gaia_dr2.ra, gaia_dr2.dec)


# Comment:
# External table from gaia. Stars for simulation and calibration.

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
    id = Column(BigInteger, nullable=False)  # Unique identifier"),
    # Fields
    ra = Column(Float(53), nullable=False)  # Right Ascension of the object (deg)"),
    dec = Column(Float(53), nullable=False)  # Declination of the object (deg)"),
    ra_err = Column(Float(24), nullable=False)  # Error in RA (arcsec)"),
    dec_err = Column(Float(24), nullable=False)  # Error in DEC (arcsec)"),
    mag_r = Column(Float(24), nullable=False)  # Model magnitude in r filter"),
    mag_b = Column(Float(24), nullable=False)  # Model magnitude in b filter"),
    mag_err_r = Column(Float(24), nullable=False)  # Mean error on mag_r"),
    mag_err_b = Column(Float(24), nullable=False)  # Mean error on mag_b"),


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
    objno = Column(BigInteger, nullable=False)  # DEEP2 object number"),
    # Fields
    ra = Column(Float(53), nullable=False)  # Right Ascension (in decimal degrees, J2000)"),
    dec = Column(Float(53), nullable=False)  # Declination (in decimal degrees, J2000),
    magb = Column(Float(24), nullable=False)  # CFHT B-band magnitude (AB) from Coil et al. 2004"),
    magr = Column(Float(24), nullable=False)  # CFHT R-band magnitude (AB) from Coil et al. 2004"),
    magi = Column(Float(24), nullable=False)  # CFHT I-band magnitude (AB) from Coil et al. 2004"),
    magberr = Column(Float(24), nullable=False)  # B-band magnitude error"),
    magrerr = Column(Float(24), nullable=False)  # R-band magnitude error"),
    magierr = Column(Float(24), nullable=False)  # I-band magnitude error"),
    rg = Column(Float(24),
                nullable=False)  # estimated R-band radius of object (sigma of Guassian fit in units of pixels --- 0.207”/pix)"),
    e2 = Column(Float(24), nullable=False)  # ellipticity defined as E2 = (1 - b/a)"),
    pa = Column(Float(24), nullable=False)  # object PA (degrees E of N)"),
    pgal = Column(Float(24),
                  nullable=False)  # the probability (0 - 1) that the sources is a galaxy for unresolved galaxies, 3 if resolved"),
    sfd_ebv = Column(Float(24), nullable=False)  # E(B-V) from Schlegel, Finkbeiner, and Davis dust map"),
    m_b = Column(Float(24), nullable=False)  # absolute B-band magnitude (AB, h = 1) from Willmer et al. (2006)"),
    ub = Column(Float(24), nullable=False)  # rest-frame U-B color (AB) from Willmer et al. (2006)"),
    objname = Column(String(8), nullable=False)  # the 8-digit DEEP2 object id (not always the same as OBJNO)"),
    mask = Column(BigInteger, nullable=False)  # the DEEP2/DEEP3 slitmask number on which the object was observed"),
    slit = Column(BigInteger, nullable=False)  # the slitlet number (on mask MASKNAME) in which the object was placed"),
    date = Column(Date, nullable=False)  # Date on which the mask was observed (YYYY-MM-DD)"),
    mjd = Column(Float(24), nullable=False)  # Modified Julian date of observation"),
    slitra = Column(Float(24), nullable=False)  # RA of slit center"),
    slitdec = Column(Float(24), nullable=False)  # Dec of slit center"),
    slitpa = Column(Float(24), nullable=False)  # PA (degrees E of N) of slit"),
    slitlen = Column(Float(24), nullable=False)  # length of slit (arcsec)"),
    z = Column(Float(24), nullable=False)  # observed best-fitting redshift"),
    zbest = Column(Float(24), nullable=False)  # best redshift (corrected for heliocentric motion)"),
    zerr = Column(Float(24), nullable=False)  # redshift error (zerr < 0 indicates problematic z fit)"),
    zquality = Column(Integer, nullable=False)  # redshift quality code, Q"),
    obj_type = Column(String(6), nullable=True)  # type of best-fitting template (e.g., GALAXY or STAR)"),
    star_type = Column(String(6), nullable=True)  # coarse classification for stellar templates"),
    rchi2 = Column(Float(24), nullable=False)  # reduced chi-squared value for the redshift fit"),
    dof = Column(BigInteger, nullable=False)  # degrees of freedom for redshift fit"),
    vdisp = Column(Float(24), nullable=False)  # velocity dispersion in km/s"),
    vdisperr = Column(Float(24), nullable=False)  # error in velocity dispersion"),
    comment = Column(String(47), nullable=True)  # comment field"),
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
    paudm_id = Column(BigInteger, nullable=False)  # PAUdm numeric id for reference table"),
    # Fields
    alpha_j2000 = Column(Float(53), nullable=False)  # Right Ascension (in decimal degrees, J2000),
    delta_j2000 = Column(Float(53), nullable=False)  # Declination (in decimal degrees, J2000),
    flag = Column(Integer, nullable=False)  # whether the object has flags or not, 0: no flags, 1: some flags,
    class_star = Column(Float(24),
                        nullable=False)  # CLASS_STAR usually gives a cleaner star sample than star_flag, but can lead to serious incompleteness in a galaxy sample,
    e1 = Column(Float(24),
                nullable=False)  # expectation values of galaxy ellipticity, from the marginalised galaxy ellipticity likelihood surface, to be used for shear measurement,
    e2 = Column(Float(24),
                nullable=False)  # expectation values of galaxy ellipticity, from the marginalised galaxy ellipticity likelihood surface, to be used for shear measurement,
    scalelength = Column(Float(24),
                         nullable=False)  # lensfit galaxy model scalelength, marginalised over other parameters, as defined by Miller et al (2012),
    bulge_fraction = Column(Float(24),
                            nullable=False)  # ratio of the ﬂux in the bulge component to the total ﬂux (often written B/T)),
    snratio = Column(Float(24),
                     nullable=False)  # signal-to-noise ratio of the object, measured on a stack of the supplied exposures within a limiting isophote 2sigma above the noise,
    z_b = Column(Float(24), nullable=False)  # most likely redshift from BPZ,
    z_b_min = Column(Float(24), nullable=False)  # z_min 95% confidence interval from BPZ,
    z_b_max = Column(Float(24), nullable=False)  # z_max 95% confidence interval from BPZ,
    t_b = Column(Float(24),
                 nullable=False)  # BPZ spectral type. 1=CWW-Ell, 2=CWW-Sbc, 3=CWW-Scd, 4=CWW-Im, 5=KIN-SB3, 6=KIN-SB2. Note that we use a recalibrated template set described in Capak et al. (2004) and that the templates are interpolated, hence fractional types occur,
    odds = Column(Float(24), nullable=False)  # likelihood that z_b is correct from BPZ,
    z_ml = Column(Float(24), nullable=False)  # maximum likelihood result (with "flat" unphysical prior) from BPZ,
    t_ml = Column(Float(24), nullable=False)  # maximum likelihood result from BPZ,
    chi_squared_bpz = Column(Float(24), nullable=False)  # "modified" chi square,
    star_flag = Column(Integer,
                       nullable=False)  # star_flag is optimized for galaxy studies, to keep an almost 100% complete galaxy sample with low (but not vanishing) stellar contamination,
    mag_u = Column(Float(24), nullable=False)  # observed magnitude in the u-band,
    magerr_u = Column(Float(24), nullable=False)  # error in the observed magnitude in the u-band,
    mag_g = Column(Float(24), nullable=False)  # observed magnitude in the g-band,
    magerr_g = Column(Float(24), nullable=False)  # error in the observed magnitude in the g-band,
    mag_r = Column(Float(24), nullable=False)  # observed magnitude in the r-band,
    magerr_r = Column(Float(24), nullable=False)  # error in the observed magnitude in the r-band,
    mag_i = Column(Float(24), nullable=False)  # observed magnitude in the i-band,
    magerr_i = Column(Float(24), nullable=False)  # error in the observed magnitude in the i-band,
    mag_y = Column(Float(24), nullable=False)  # observed magnitude in the y-band,
    magerr_y = Column(Float(24), nullable=False)  # error in the observed magnitude in the y-band,
    mag_z = Column(Float(24), nullable=False)  # observed magnitude in the z-band,
    magerr_z = Column(Float(24), nullable=False)  # error in the observed magnitude in the z-band,


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
    mod_gal = Column(Integer,
                     nullable=False)  # best fit template  (1->8 Ell-S0, 9->15 Sa-Sc, 16->19 Sd-Sdm from Polletta et al., >=20 from BC03).
    type = Column(Integer, nullable=False)  # 0 if gal, 1 if star, 2 if Xray, 3 if IA>25.5, -99 if masked
    Umag = Column(Float(24), nullable=False)  # CFHT u*   aperture magnitude (3" aperture)        3911-0 +- 538/2
    Bmag = Column(Float(24), nullable=False)  # Subaru B   aperture magnitude (3" aperture)       4439.6 +- 806.7/2
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
    acs_a_image = Column(Float(24), nullable=False)  # additional extinction  E(B-V)
    acs_b_image = Column(Float(24), nullable=False)  # additional extinction  E(B-V)
    acs_theta_image = Column(Float(24), nullable=False)  # additional extinction  E(B-V)
    acs_mag_auto = Column(Float(24), nullable=False)  # additional extinction  E(B-V)
    acs_magerr_auto = Column(Float(24), nullable=False)  # additional extinction  E(B-V)
    zspec = Column(Float(24), nullable=False)  # spectroscopic redshift from ADP.2015-10-15T12_55_10.777_t2.txt
    conf = Column(Float(24), nullable=False)  # confidence class from ADP.2015-10-15T12_55_10.777_t2.txt
    F814W = Column(Float(24), nullable=False)  # HST  target magnitude from ADP.2015-10-15T12_55_10.777_t2.txt
    zfits = Column(String(60), nullable=True)  # filename of associated 1-d spectra ADP.2015-10-15T12_55_10.777_t2.txt
    r50 = Column(Float(24), nullable=False)  # ZEST semi-major axis length of ellipse encompassing 50% of total light
    sersic_n_gim2d = Column(Float(24), nullable=False)  # GIM2D Sersic index


Index('ik_cosmoslocation', COSMOS.ra, COSMOS.dec)


# Comment:
# External table from zCOSMOS (DR3). Sources with accurate redshifts for forced photometry and validation.


# Spectra-Convolved Table
class SpecConv(Base):
    __tablename__ = 'spec_conv'
    __table_args__ = (
        # Constraints
        PrimaryKeyConstraint('id'),
        UniqueConstraint('spec_id', 'spec_cat', 'band', 'instrument'),
    )
    # Keys
    id = Column(BigInteger, nullable=False)  # Unique identifier"),
    # Fields
    spec_id = Column(BigInteger, nullable=False)  # Unique identifier for spectra source
    spec_cat = Column(String(16), nullable=True)  # Spectra catalogue name [COSMOS, SDSS, DEEP2,...]
    band = Column(String(16), nullable=True)  # Band name to be convolved
    instrument = Column(String(16), nullable=True)  # Instrument of the band name [PAU, SDSS, DES,...]
    # fluxes
    flux = Column(Float(24), nullable=True)  # fluxed measured in the spectra convolved with the specified filter
    flux_err = Column(Float(24), nullable=True)  # associated flux error

Index('ik_specconv_refid', SpecConv.spec_cat, SpecConv.spec_id)

# Comment:
# Contains the convolved fluxes derived from spectra observations from external surveys (i.e. SDSS, COSMOS, DEEP2...)


# Match to Spectra table
class MatchToSpec(Base):
    __tablename__ = 'match_to_spec'
    __table_args__ = (
        # Constraints
        PrimaryKeyConstraint('ref_id', 'ref_cat', 'spec_id', 'spec_cat'),
    )
    # Keys
    ref_id = Column(BigInteger, nullable=False)  # Unique identifier for reference catalogue
    ref_cat = Column(String(16), nullable=False)  # Reference catalogue name
    spec_id = Column(BigInteger, nullable=False)  # Unique identifier for spectra reference catalogue
    spec_cat = Column(String(16), nullable=False)  # Reference spectra catalogue name


Index('ik_matchtospec_ref', MatchToSpec.ref_cat, MatchToSpec.ref_id)


### SIMULATION TABLES ###

# Truth Objects table
class Truth_Object(Base):
    __tablename__ = 'truth_object'
    __table_args__ = (
        # Constraints
        PrimaryKeyConstraint('id'),
        ForeignKeyConstraint(['production_id'], ['production.id'], onupdate='CASCADE', ondelete='CASCADE'),
    )
    # Keys
    id = Column(BigInteger, nullable=False)  # Unique identifier"),
    # Fields
    production_id = Column(Integer, nullable=False)  # Production number"),
    stargalaxy = Column(Boolean, nullable=False)  # 0(False):galaxy, 1(True):star"),
    sed_type = Column(Float(24), nullable=False)  # SED type"),
    sed_em_line = Column(Float(24), nullable=False)  # SED emission line type"),
    ra = Column(Float(53), nullable=False)  # Truth Right AscensSion position"),
    dec = Column(Float(53), nullable=False)  # Truth Declination position"),
    mag_u = Column(Float(24), nullable=False)  # Truth simulated magnitude for filter u"),
    mag_g = Column(Float(24), nullable=False)  # Truth simulated magnitude for filter g"),
    mag_r = Column(Float(24), nullable=False)  # Truth simulated magnitude for filter r"),
    mag_i = Column(Float(24), nullable=False)  # Truth simulated magnitude for filter i"),
    mag_z = Column(Float(24), nullable=False)  # Truth simulated magnitude for filter z"),
    mag_y = Column(Float(24), nullable=False)  # Truth simulated magnitude for filter y"),
    mag_NB455 = Column(Float(24), nullable=False)  # Truth simulated magnitude for PAU narrow filter 01"),
    mag_NB465 = Column(Float(24), nullable=False)  # Truth simulated magnitude for PAU narrow filter 02"),
    mag_NB475 = Column(Float(24), nullable=False)  # Truth simulated magnitude for PAU narrow filter 03"),
    mag_NB485 = Column(Float(24), nullable=False)  # Truth simulated magnitude for PAU narrow filter 04"),
    mag_NB495 = Column(Float(24), nullable=False)  # Truth simulated magnitude for PAU narrow filter 05"),
    mag_NB505 = Column(Float(24), nullable=False)  # Truth simulated magnitude for PAU narrow filter 06"),
    mag_NB515 = Column(Float(24), nullable=False)  # Truth simulated magnitude for PAU narrow filter 07"),
    mag_NB525 = Column(Float(24), nullable=False)  # Truth simulated magnitude for PAU narrow filter 08"),
    mag_NB535 = Column(Float(24), nullable=False)  # Truth simulated magnitude for PAU narrow filter 09"),
    mag_NB545 = Column(Float(24), nullable=False)  # Truth simulated magnitude for PAU narrow filter 10"),
    mag_NB555 = Column(Float(24), nullable=False)  # Truth simulated magnitude for PAU narrow filter 11"),
    mag_NB565 = Column(Float(24), nullable=False)  # Truth simulated magnitude for PAU narrow filter 12"),
    mag_NB575 = Column(Float(24), nullable=False)  # Truth simulated magnitude for PAU narrow filter 13"),
    mag_NB585 = Column(Float(24), nullable=False)  # Truth simulated magnitude for PAU narrow filter 14"),
    mag_NB595 = Column(Float(24), nullable=False)  # Truth simulated magnitude for PAU narrow filter 15"),
    mag_NB605 = Column(Float(24), nullable=False)  # Truth simulated magnitude for PAU narrow filter 16"),
    mag_NB615 = Column(Float(24), nullable=False)  # Truth simulated magnitude for PAU narrow filter 17"),
    mag_NB625 = Column(Float(24), nullable=False)  # Truth simulated magnitude for PAU narrow filter 18"),
    mag_NB635 = Column(Float(24), nullable=False)  # Truth simulated magnitude for PAU narrow filter 19"),
    mag_NB645 = Column(Float(24), nullable=False)  # Truth simulated magnitude for PAU narrow filter 20"),
    mag_NB655 = Column(Float(24), nullable=False)  # Truth simulated magnitude for PAU narrow filter 21"),
    mag_NB665 = Column(Float(24), nullable=False)  # Truth simulated magnitude for PAU narrow filter 22"),
    mag_NB675 = Column(Float(24), nullable=False)  # Truth simulated magnitude for PAU narrow filter 23"),
    mag_NB685 = Column(Float(24), nullable=False)  # Truth simulated magnitude for PAU narrow filter 24"),
    mag_NB695 = Column(Float(24), nullable=False)  # Truth simulated magnitude for PAU narrow filter 25"),
    mag_NB705 = Column(Float(24), nullable=False)  # Truth simulated magnitude for PAU narrow filter 26"),
    mag_NB715 = Column(Float(24), nullable=False)  # Truth simulated magnitude for PAU narrow filter 27"),
    mag_NB725 = Column(Float(24), nullable=False)  # Truth simulated magnitude for PAU narrow filter 28"),
    mag_NB735 = Column(Float(24), nullable=False)  # Truth simulated magnitude for PAU narrow filter 29"),
    mag_NB745 = Column(Float(24), nullable=False)  # Truth simulated magnitude for PAU narrow filter 30"),
    mag_NB755 = Column(Float(24), nullable=False)  # Truth simulated magnitude for PAU narrow filter 31"),
    mag_NB765 = Column(Float(24), nullable=False)  # Truth simulated magnitude for PAU narrow filter 32"),
    mag_NB775 = Column(Float(24), nullable=False)  # Truth simulated magnitude for PAU narrow filter 33"),
    mag_NB785 = Column(Float(24), nullable=False)  # Truth simulated magnitude for PAU narrow filter 34"),
    mag_NB795 = Column(Float(24), nullable=False)  # Truth simulated magnitude for PAU narrow filter 35"),
    mag_NB805 = Column(Float(24), nullable=False)  # Truth simulated magnitude for PAU narrow filter 36"),
    mag_NB815 = Column(Float(24), nullable=False)  # Truth simulated magnitude for PAU narrow filter 37"),
    mag_NB825 = Column(Float(24), nullable=False)  # Truth simulated magnitude for PAU narrow filter 38"),
    mag_NB835 = Column(Float(24), nullable=False)  # Truth simulated magnitude for PAU narrow filter 39"),
    mag_NB845 = Column(Float(24), nullable=False)  # Truth simulated magnitude for PAU narrow filter 40"),
    mag_n41 = Column(Float(24), nullable=True)  # Truth simulated magnitude for PAU narrow filter 41"),
    mag_n42 = Column(Float(24), nullable=True)  # Truth simulated magnitude for PAU narrow filter 42"),
    bulge2flux_ratio = Column(Float(24), nullable=True)  # Galaxy only - Bulge to total flux ratio"),
    bulge_length = Column(Float(24), nullable=True)  # Galaxy only - Bulge length scale (arcsec)"),
    bulge_axis_ratio = Column(Float(24), nullable=True)  # Galaxy only - Bulge projected axis ratio"),
    bulge_angle = Column(Float(24), nullable=True)  # Galaxy only - Bulge position angle"),
    disk_length = Column(Float(24), nullable=True)  # Galaxy only - Disk length scale (arcsec)"),
    disk_axis_ratio = Column(Float(24), nullable=True)  # Galaxy only - Disk projected axis ratio"),
    disk_angle = Column(Float(24), nullable=True)  # Galaxy only - Disk position angle"),
    redshift = Column(Float(24), nullable=True)  # Galaxy only - Truth redshift"),
    # Constraints
    production = relationship('Production', back_populates="truth_objects")


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
    id = Column(Integer, nullable=False)  # Unique identifier"),
    production_id = Column(Integer, nullable=False)  # Production number"),
    # Fields
    exp_num = Column(Integer, nullable=True)  # comment="Camera exposure unique identifier"),
    field = Column(String(16), nullable=True)  # comment="Field Name"),
    ra = Column(Float(53), nullable=True)  # Target RA"),
    dec = Column(Float(53), nullable=True)  # Target DEC"),
    filtertray = Column(String(16), nullable=True)  # Filter Tray name"),
    kind = Column(Enum('BIAS', 'FLAT', 'TARGET', name='target_kind'), nullable=False)  # Target kind"),
    status = Column(Enum('PLANNED', 'SCHEDULED', 'SIMULATED', 'EXPOSED', name='target_status'),
                    nullable=False)  # comment="Target status"),
    rjd_obs = Column(Float(53), nullable=True)  # Observation Reduced Modified Julian Day"),
    exptime = Column(Float(24), nullable=True)  # comment="Exposure Time",
    airmass = Column(Float(24), nullable=True)  # comment="Airmass",
    seeing = Column(Float(24), nullable=True)  # comment="Seeing",
    moon_mag = Column(Float(24), nullable=True)  # comment="Moon Magnitude",
    moon_phase = Column(Float(24), nullable=True)  # comment="Moon Phase",
    moon_distance = Column(Float(24), nullable=True)  # comment="Moon Distance (degrees)",
    moon_set = Column(Boolean, nullable=True)  # comment="Moon is set?",
    ecliptic_dist = Column(Float(24), nullable=True)  # comment="Distance to ecliptic (degrees)",
    ecliptic_zodiacal = Column(Float(24), nullable=True)  # comment="Ecliptic Zodiacal Light",
    ecliptic_airglow = Column(Float(24), nullable=True)  # comment="Ecliptic Airglow Light",

    # Relationships
    bkg_mags = relationship('Bkg_mag', back_populates="target")
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
    target_id = Column(Integer, nullable=False)  # comment="Unique identifier",
    filter = Column(String(16), nullable=False)  # comment="Filter Name"),
    # Fields
    mag = Column(Float(24), nullable=False)  # comment="Magnitude value"),
    # Documentation
    # comment="Background Magnitude",

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
        ForeignKeyConstraint(['job_id'], ['job.id'], ondelete='CASCADE')

    )
    # Keys
    id = Column(Integer, nullable=False)  # Unique identifier"),
    # Fields
    job_id = Column(Integer, nullable=False)  # Job identifier"),
    ref = Column(String(32), nullable=False)  # QC Reference code"),
    check_name = Column(String(64), nullable=False)  # QC check name"),
    min_value = Column(Float(24), nullable=False)  # Minimum Range value"),
    max_value = Column(Float(24), nullable=False)  # Maximum Range value"),
    value = Column(Float(24), nullable=False)  # Measured value"),
    units = Column(String(24), nullable=False)  # Measurement units"),
    qc_pass = Column(Boolean, nullable=False)  # QC pass?"),
    time = Column(DateTime, nullable=False, default=func.current_timestamp())  # comment="Timestamp of creation", ),
    plot_file = Column(Text, nullable=True)  # comment="Plot file full name and path"),

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

    for table_name, columns in criteria.items():
        for column_name, filtr in columns.items():
            if isinstance(filtr, list):
                q = q.filter(tables[table_name].columns[column_name].in_(filtr))
            elif isinstance(filtr, dict):
                for operator, value in filtr.items():
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


def main(argv=None):
    if not argv:
        argv = sys.argv

    (options,) = opt_parser(argv)[:1]

    init(options.url)

    if not options.force:
        answer = input("The database will be ERASED and recreated. Do you want to proceed? (y/N): ")
        if answer.upper() != "Y":
            return

    recreate()


def recreate():
    # Drop and recreate the database
    metadata.drop_all()
    # metadata.create_all()
    if metadata.bind.url.database == 'paudb' and metadata.bind.url.drivername == 'postgresql':
        metadata.create_comments()
        metadata.set_permissions()


if __name__ == '__main__':
    import sys

    sys.exit(main(sys.argv))
