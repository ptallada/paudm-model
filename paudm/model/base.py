#! /usr/bin/env python
# -*- coding: utf-8 -*-

import sqlalchemy

from sqlalchemy.orm import sessionmaker

class Table(sqlalchemy.Table):
    def __new__(cls, *args, **kwargs):
        # Remove comment to avoid exception
        kwargs.pop('comment', None)
        return super(Table, cls).__new__(cls, *args, **kwargs)
    
    def __init__(self, *args, **kwargs):
        doc = kwargs.pop('comment', None)
        super(Table, self).__init__(*args, **kwargs)
        self.__doc__ = doc

class Column(sqlalchemy.Column):
    def __init__(self, *args, **kwargs):
        doc = kwargs.pop('comment', None)
        super(Column, self).__init__(*args, **kwargs)
        self.__doc__ = doc

class MetaData(sqlalchemy.MetaData):
    def create_comments(self, bind=None):
        if bind is None:
            bind = sqlalchemy.schema._bind_or_error(self)
        session = sessionmaker(bind)()
        for t in self.sorted_tables:
            if t.__doc__:
                session.execute("COMMENT ON TABLE \"%s\" IS '%s'"
                    % (t.name, t.__doc__.replace("'", "''").strip()))
            for c in t.columns:
                if c.__doc__:
                    session.execute("COMMENT ON COLUMN \"%s\".\"%s\" IS '%s'"
                        % (t.name, c.name, c.__doc__.replace("'", "''").strip()))
        # session.commit()
        
    def set_permissions(self, bind=None):
        if bind is None:
            bind = sqlalchemy.schema._bind_or_error(self)
        session = sessionmaker(bind)()
        for t in self.sorted_tables:
            session.execute('ALTER TABLE "%s" OWNER TO admin' % t.name)
            session.execute('GRANT ALL ON TABLE "%s" TO postgres' % t.name)
            session.execute('GRANT ALL ON TABLE "%s" TO admin' % t.name)
            session.execute('GRANT SELECT, UPDATE, INSERT ON TABLE "%s" TO production' % t.name)
            session.execute('GRANT SELECT, UPDATE, INSERT ON TABLE "%s" TO analysis' % t.name)
        for t in ['mosaic_kind', 'target_kind', 'event_status', 'job_status', 'target_status']:
            session.execute('ALTER TYPE "%s" OWNER TO admin' % t)
        for t in ['event', 'job', 'job_dependency']:
            session.execute('GRANT SELECT, UPDATE, INSERT ON TABLE "%s" TO events' % t)
        # session.commit()
