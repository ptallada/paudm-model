"""
Collection of function to operate the downoload and the load of the files from and to the storage
to and from the working directory.

The environment defined in the common configuration defines which commands to use for the copy:

PC environment: define function here
GRID environm.: call to tools/grid functions for download


Environment should be pass as argument
Here only general functions
"""

import os
import shutil

import logging
log = logging.getLogger('paudm.tools.common.archive')

from paudm.tools.grid import storage_utils


# Makedirs wrapper to avoid race conditions
def mkdir_p(path):
    try:
        os.makedirs(path)
    except os.error, e:
        if e.errno != errno.EEXIST:
            raise

def list_directory(path, env_mode, basename):
    """
    Make list of dictionaries for the FILES in the specified path.
    Dictionaries keywords: name, path
    """
    # config['instrument']['BASENAME'] needed
    # config['general']['ENV_MODE']  needed
    # List files
    file_list = []
    if env_mode != 'PC':
        # In storage
        file_list = storage_utils.list_storage_directory(path = path, pattern = basename)
    else:
        # In local path
        dir_list = os.listdir(path)
        for filename in dir_list:
            file_list.append({'name': filename, 'path': path})
    
    # No files found
    if len(file_list) == 0:
        log.warning("No files found in storage path.")
        
    return file_list


def list_query(query):
    """
    Make list of file type dictionaries for FILES in query.
    Dictionaries keywords: name, path
    """
    
    file_list = []
    
    # Create dictionary of mosaic files
    for db_element in query:
        file_list.append({'name': db_element.filename, 'path': db_element.archivepath})
        
    # No mosaics found
    if len(file_list) == 0:
        log.warn("No files found for current query.")
    
    return file_list
    

def download(config_grid, file_list, dest_path, env_mode , grid_protocol = None):
    #config['general']['ENV_MODE'] needed
    #config['grid']['PROTOCOL'] optional
    # Add only files not yet in the working directory
    file_list_to_download = []
    for file in file_list:
        if os.path.isdir(os.path.join(dest_path, file['name'])):
            log.warning("Filename is a directory")
            continue
        if not os.path.exists(os.path.join(dest_path, file['name'])):
            log.debug("file %s not in working path. Scheduled to download."% (file['name']))
            file_list_to_download.append(file)
    
    # PRESTAGE
    if env_mode != 'PC':
        # If environment is GRID, download from storage
        grid_protocol = "SRM" if grid_protocol==None else grid_protocol
        storage_utils.download(config_grid, dest_path = dest_path, file_list = file_list_to_download, protocol = grid_protocol)
    else:
        # If environment is PC, just do symbolic link
        for file in file_list_to_download:
            os.symlink(os.path.join(file['path'], file['name']), os.path.join(dest_path, file['name']))
    log.info("Download to input path completed")
    


def copy_to_nfs(file_in, file_out):
    
    if file_in == None or not os.path.exists(file_in):
        msg = "No file %s to copy to NFS"% file_in
        log.error(msg)
        raise Exception, msg

    if file_out == None:
        msg = "copy_to_nfs: No NFS destination given"
        log.error(msg)
        raise Exception, msg
    
    if not os.path.exists(os.path.dirname(file_out)):
        directory = os.path.dirname(file_out)
        log.info( "Making directory %s" %directory )
        mkdir_p(directory)
        os.chmod(directory, 0775)
            
    log.info("copy %s to %s"% (file_in,file_out))
    if os.path.exists(file_out):
        os.remove(file_out)
    shutil.copyfile(file_in, file_out)
    


def upload(config, file_in = None, file_out = None):
    
    if file_in == None or not os.path.exists(file_in):
        log.warning("No file to upload... file = %s"% file_in)
        return
    
    # PC mode
    if config['general']['env_mode'] == 'PC':
        # Make whatever directories are necessary.  Archive path should exist.
        directory = os.path.dirname(file_out)
        if not os.path.exists(directory):
            log.info( "Making directory %s" %directory )
            mkdir_p(directory)
            os.chmod(directory, 0775)
        
        log.info("Copying file %s to archive at %s..." %(file_in, file_out))
        shutil.move(file_in, file_out)
        log.info("file copy completed successfully.")
    
    # GRID mode
    else:
        archive_type = config['common']['grid']['archive_type']
        if archive_type == 'STORAGE':
            from paudm.tools.grid.storage_utils import upload as grid_upload
            log.info("Uploading file %s to archive in STORAGE at location %s..." %(file_in, file_out))
            grid_upload(config['common'], file_in = file_in, file_out = file_out)
            log.info("file upload completed successfully.")
        
        elif archive_type == 'NFS':
            filename = os.path.basename(file_in)
            log.info("Copying file %s to NFS..." %(file_in))
            if file_out == None:
                file_out = os.path.join(config['common']['grid']['NFS_PATH'], 'sandbox_data', filename)
            copy_to_nfs(file_in, file_out)
            log.info("file copy completed successfully.")
            
        elif archive_type == 'NONE':
            log.info("Archive type is NONE. Will not upload the file: %s" %file_in)
            return 0
        else:
            msg = "Incorrect archive type %s for GRID. Set ARCHIVE, NFS or NONE." %archive_type
            log.error(msg)
            raise Exception, msg
            
    log.info("Upload to destination path completed")
    


