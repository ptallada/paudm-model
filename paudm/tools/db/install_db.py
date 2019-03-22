#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This script is intended to be run to install a new database

import sys
import argparse


def _parse_args(args):
    parser = argparse.ArgumentParser(prog='install_db', argument_default=argparse.SUPPRESS,
                                     epilog="To recreate the complete db and eliminate the private SDSS and USNO tables: python install_db.py -D -d SDSS USNO")

    parser.add_argument('-D', '--drop-all-before-create', default=False, dest='drop_all', action='store_true',
                        help='Create a brand new database')

    parser.add_argument('-d', '--drop-tables', default=None, metavar='TABLENAME', nargs='+',
                        help='The list of private tables that will be deleted, in order to access the public ones.')

    options = vars(parser.parse_args(args))
    return options


def main(argv=None):
    from paudm.tools.db import model
    from paudm.environment import context

    print("Welcome to the install_db script.")
    if not argv:
        argv = sys.argv[1:]

    options = _parse_args(argv)
    drop_all = options.pop('drop_all')
    drop_tables = options.pop('drop_tables')

    with context.init(None, ""):
        session = context.get_session()

        # Create DB
        # Initializing context
        session = None
        if drop_all:
            answer = input("The database will be ERASED and recreated. Do you want to proceed? (y/N): ")
            if answer.upper() == "Y":
                model.metadata.drop_all()
                print("Database successfully erased")
            else:
                print("Database not erased.")

        model.metadata.create_all()

        if drop_tables:
            answer = input(
                "The private tables %s will be deleted, in order to access the public catalogs. Do you want to proceed? (y/N): " % str(
                    drop_tables))
            if answer.upper() == "Y":
                for table in drop_tables:
                    model.metadata.tables[table.lower()].drop()
                    print("Private table %s deleted. You will be able to access the public one" % str(table))
            else:
                print("Your private tables will be preserved and  you will be able to write on them.")

    print("...DB install finishes here.")


if __name__ == '__main__':
    sys.exit(main())
