# -*- coding: utf-8 -*-
"""
Created on Wed Apr 24 08:05:25 2013

@author: dpiscia & serrano
"""
import datetime
from paudm.tools.db import model
from sqlalchemy import or_
from sqlalchemy import func


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
        log.warn("Quality Control %s not passed: %f %s" % (qc_ref, value, qc_constants.qc[qc_ref]['Units']))


def update_tree(session, job, qc_pass):
    if qc_pass == "":
        quality_control_value = session.query(model.Quality_control).filter(
            model.Quality_control.job_id == job.id,
            model.Quality_control.check_name != "general",
            model.Quality_control.qc_pass == False,
        ).first() is not None

        qc = model.session.query(model.Quality_control).filter(
            model.Quality_control.check_name == "general",
            model.Quality_control.job_id == job.id,
        ).first()

        if qc:
            # Should it be deleted, or there can be some case where it exits and it is correct??
            session.delete(qc)

        else:
            qc = model.Quality_control(
                job_id=job.id,
                ref="general",
                check_name="general",
                min_value=0,
                max_value=0,
                value=0,
                units=0,
                qc_pass=quality_control_value,
                time=datetime.datetime.now(),
                plot_file="",
            )
            session.add(qc)
            update_tree(session, job, quality_control_value)
    else:
        qc = session.query(model.Quality_control).filter(
            model.Quality_control.check_name == "general",
            model.Quality_control.job_id == job.id
        ).first()

        if qc:
            qc.qc_pass = (qc.qc_pass and qc_pass)

        else:
            qc = model.Quality_control(
                job_id=job.id,
                ref="general",
                check_name="general",
                min_value=0,
                max_value=0,
                value=0,
                units=0,
                qc_pass=qc_pass,
                time=datetime.datetime.now(),
                plot_file=""
            )

            session.add(qc)

        if job.superjob:
            update_tree(session, job.superjob, qc_pass)


def update_parent_qc(session, job):
    if not job.superjob:
        combined_qc_pass = session.query(func.bool_and(model.Quality_control.qc_pass))\
            .filter(model.Quality_control.job_id == job.id).one()
    else:
        combined_qc_pass = session.query(func.bool_and(model.Quality_control.qc_pass))\
            .filter(or_(model.Quality_control.job_id == job.id, model.Quality_control.job_id == job.superjob.id)).one()

    qc = model.Quality_control(
        job_id=job.id,
        ref="general",
        check_name="general",
        min_value=0,
        max_value=0,
        value=0,
        units=0,
        qc_pass=combined_qc_pass[0],
        time=datetime.datetime.now(),
        plot_file="",
    )
    session.add(qc)




