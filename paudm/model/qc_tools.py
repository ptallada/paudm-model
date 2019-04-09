# -*- coding: utf-8 -*-
"""
Created on Wed Apr 24 08:05:25 2013

@author: dpiscia & serrano
"""
import datetime
from . import model

import logging
log = logging.getLogger('paudm.tools.common.qc_tools')


def quality_control_entry(job_id, qc_ref, value, qc_constants, plot_name='', label=''):
    qc_pass = bool(qc_constants.qc[qc_ref]['Min_value'] <= value <= qc_constants.qc[qc_ref]['Max_value'])
    qc_entry = model.Quality_control(job_id=job_id,
                                     ref=qc_ref + label,
                                     check_name=qc_constants.qc[qc_ref]['Name'],
                                     min_value=qc_constants.qc[qc_ref]['Min_value'],
                                     max_value=qc_constants.qc[qc_ref]['Max_value'],
                                     value=value,
                                     units=qc_constants.qc[qc_ref]['Units'],
                                     qc_pass=qc_pass,
                                     time=datetime.datetime.now(),
                                     plot_file=plot_name)
    model.session.add(qc_entry)
    if not qc_pass:
        log.warning("Quality Control %s not passed: %f %s" % (qc_ref, value, qc_constants.qc[qc_ref]['Units']))


def update_parent_qc(session, job):

    combined_qc_pass = session.execute('select bool_and(qc.qc_pass) from quality_control as qc '
                                       'join job on job.id=qc.job_id '
                                       'where job.id =%d or job.super_id=%d' % (job.id, job.id)).fetchall()[0][0]

    if combined_qc_pass is not None:
        qc = model.Quality_control(
            job_id=job.id,
            ref="general",
            check_name="general",
            min_value=0,
            max_value=0,
            value=0,
            units=0,
            qc_pass=combined_qc_pass,
            time=datetime.datetime.now(),
            plot_file="",
        )
        session.add(qc)




