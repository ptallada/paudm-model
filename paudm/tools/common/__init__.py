# COMMON initialization

import os
import sys
import time
import errno
import logging.config
import yaml
log = logging.getLogger('paudm.tools.common')

from paudm.tools.db import model
from pkg_resources import resource_string
_config = None

#print "Initializing PAUdm system..."

def get_config():
 
    global _config

    if _config != None:
        return _config
    
    ### Load Common Configuration ###
    print "Loading common configuration..."
    config = yaml.safe_load(resource_string('paudm.tools.common','config/common.yaml'))
    
    ### Load Command Line configuration (overwrites default) ###
    print "Loading command line configuration..."
    
    # Get executed python file
    python_file = ''
    if len(sys.argv) > 0:
        python_file = sys.argv.pop(0)
    # Only get pipeline name when calling shuttle
    if python_file.rfind('shuttle') != -1:
        if len(sys.argv) != 0:
            # Pipeline specified
            if not (sys.argv[0].startswith('+') or sys.argv[0].startswith('-')):
                run_pipeline = sys.argv.pop(0)
                valid_pipelines = ['pixelsim', 'nightly', 'memba', 'analysis', 'common', 'commissioning']
                if run_pipeline not in valid_pipelines:
                    print >> sys.stderr, "Invalid pipeline name. The pipelines are: "
                    for pipeline in valid_pipelines:
                        print >> sys.stderr, " - " + pipeline
                    sys.exit()
                # Get run pipeline 
                config['RUN_PIPELINE'] = run_pipeline
                print "Shuttle will launch the %s pipeline!" %run_pipeline
            # Pipeline NOT specified
            else:
                config['RUN_PIPELINE'] = None
        # Pipeline NOT specified
        else:
            config['RUN_PIPELINE'] = None
            
    # Get configuration options
    options = sys.argv
    for option in options:
        if option.startswith('+'):
            try:
                option = option[1:]
                location, value_string = option.split('=')
                section, variable = location.split(',')
                value_converted = yaml.safe_load(value_string)
                print "Overwriting COMMON variable %s in %s section with value: "%(variable, section) + str(value_converted)
                config[section][variable] = value_converted
            except Exception:
                print >> sys.stderr, "Invalid command line configuration: %s" %option
    
    
    # Load Instrument Constants
    general_instr_file = open('paudm/resources/instrument/' + config['general']['PROJECT'].lower() + '/general.yaml')
    config['instrument'] = yaml.safe_load(general_instr_file)
    
    ### Configure Users MASK ###
    print "Setting file mode creation mask to 0002..."
    os.umask(0002)
    
    if config['general']['WORKING_PATH'][0] == '$':
        working_path = os.getenv(config['general']['WORKING_PATH'][1:])
        # When WORKING_PATH is the name of a variable, substitute it with its value (only if exists)
        if working_path != None:
            config['general']['WORKING_PATH'] = working_path

    _config = config
    return _config


def get_logger():
    config = get_config()
    
    ### Configure LOGGING ###
    print "Setting up logging system..."
    if config['logging']['path'].lower() != 'default' :
        # Custom log file path and name (for jobs!)
        log_filename = config['logging']['path']
    else:
        # Default log file name and path
        log_filename = os.path.join(config['general']['working_path'], 'PAUdm.log')
    if not os.path.exists(log_filename):
        open(log_filename, 'w').close()
        os.chmod(log_filename, 0664)
    logging.config.fileConfig('pipeline/common/config/log.conf', {'file_name': log_filename,
                                                                  'file_level': config['logging']['file_level'],
                                                                  'console_level': config['logging']['sceen_level'],
                                                                  'sql_level': config['logging']['sql_level']
                                                                  })
    # Get logger
    log = logging.getLogger('PAUdm_init')
    log.propagate = True
    log.info("Log file in %s" %log_filename)
    return log
    
def init_db(config):
    #config = get_config()
    #log = get_logger()
    ### DATA BASE Configuration ###
    if config['general']['env_mode'] == 'PC':
        log.info("Setting up local %s SQLite database..." %config['common']['database']['db_name'])
        db_full_path = os.path.join(config['common']['database']['path'], config['common']['database']['db_name'] + '.db')
        db_config_file = resource_string('paudm.tools.db.config','local.yaml')
    else:
        log.info("Setting up PIC's %s PostrgreSQL database..." %config['common']['database']['db_name'])
        db_full_path = config['common']['database']['db_name']
        db_config_file = resource_string('paudm.tools.db.config','paudb.yaml')
        
    # Initialize DB
    settings = yaml.safe_load(db_config_file, "r")
    url      = settings.get('engine', 'sqlite:///') + db_full_path
    model.init(url)
    # model.metadata.create_all()
    # Create DB if requested
    if config['common']['database']['create_db']:
        answer = raw_input("The database will be ERASED and recreated. Do you want to proceed? (y/N): ")
        if answer.upper() != "Y":
            log.info("Data Base not restored.")
        else:
            log.info ("Creating new DB...")
            model.recreate()
            log.info("...DB successfully created.")
    
    ### Obtain REVISION Number ###
    # Will work only under svn revision
    svnversion_string = os.popen('svnversion -n').read().split(':')[0]
    try:
        while svnversion_string[-1].isalpha():
            svnversion_string = svnversion_string[:-1]
    except Exception:
        svnversion_string = '-1'
    config['general']['red_revision'] = int(svnversion_string)
    log.info ("SVN Revision number: %d" %config['general']['red_revision'])
    
    log.info("PAUdm initialization complete.")
    return model.session
