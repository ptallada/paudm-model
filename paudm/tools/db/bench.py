#! /usr/bin/env python
# -*- coding: utf-8 -*-

import datetime
import random

import model

def recreate():
    model.metadata.drop_all()
    model.metadata.create_all()
    model.metadata.create_comments()

def bench_insert(count = 100, *args):
    count = int(count)
    
    randf = random.random()
    randi = random.randint(0, 0x7fffffffL)
    
    run = model.Run(
        instrument = 'BENCH',
        log = 'BENCH',
        rjd_start = randf,
        rjd_stop = randf + 1,
        night = datetime.date.today(),
        notes = 'BENCH',
        operator = 'BENCH',
        run_number = randi,
    )
    
    production = model.Production(
        project = 'BENCH',
        origin = 'BENCH',
        origin_release = 'BENCH',
        red_revision = randi,
        red_release = 'BENCH'
    )

    mosaic = model.Mosaic(
        production = production,
        run = run,
        filename = 'BENCH %s' % randi,
        archivepath = 'BENCH %s' % randi,
        kind = random.choice(['ARC', 'BIAS', 'DARK', 'FLASH', 'FLAT', 'FOCUS', 'GLANCE', 'PUPIL', 'SCRATCH', 'SKY', 'TARGET','MBIAS','MFLAT','RED_SCI','RED_WEIGHT','RED_MASK']),
        exp_num = randi,
        obs_title = 'BENCH',
        ra = randf,
        dec = randf + 1,
        equinox = randf - 1,
        date_obs = datetime.date.today(),
        time_obs = datetime.datetime.today().time(),
        rjd_obs = randf,
        date_creat = datetime.date.today(),
        time_creat = datetime.datetime.today().time(),
        exp_time = randf,
        air_mass = randf,
        telfocus = randf,
        instrument = 'BENCH',
        filtertray = 'BENC',
        filtertray_tmp = randf,
        nextend = randi % 32767,
        guide_enabled = True,
        guide_fwhm = randf,
        guide_var = randf,
        extinction = randf,
        extinction_err = randf,
        wind_spd = randf,
        wind_dir = randf,
        amb_temp = randf,
        humidity = randf,
        pressure = randf,
        eqa_1 = randf,
        eqa_2 = randf,
        eqa_3 = randf,
        eqa_4 = randf,
        eqa_5 = randf,
    )
    
    image = model.Image(
        mosaic = mosaic,
        image_num = randi % 32767,
        ccd_num = randi % 32767,
        amp_num = randi % 32767,
        filter = 'BENC',
        wavelength = randf,
        waveband = randf,
        gain = randf,
        rdnoise = randf,
        naxis1 = randi % 32767,
        naxis2 = randi % 32767,
        zp_nightly = randf,
        zp_nightly_sigma = randf,
        psf_fwhm = randf,
        cqa_1 = randf,
        cqa_2 = randf,
        cqa_3 = randf,
        cqa_4 = randf,
        cqa_5 = randf
    )
    
    for i in range(count):
        randf = random.random()
        randi = random.randint(0, 0xffffffffL)
        
        model.session.add(model.Detection(
            production = production,
            image = image,
            band = random.choice('grizy'),
            background = randf,
            class_star = randf,
            flux_auto = randf,
            flux_err_auto = randf,
            mag_auto = randf,
            mag_err_auto = randf,
            flags = randi % 32767,
            elongation = randf,
            dec = randf,
            ra = randf + 1,
            x = randf,
            y = float(i),
            zp_offset = randf
        ))
    
    # model.session.commit()

if __name__ == '__main__':
    import sys
    model.init(env = 'BENCH')
    
    # Run bench
    globals()[sys.argv[1]](*sys.argv[2:])