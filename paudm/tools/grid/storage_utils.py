'''
Functions to be executed when running in the grid environment.

Files at PIC are permanently stored in tape. In order to access them, they should be copied in the /nfs.

Steps to perform a prestage:
 1. Locate where the original files are and where they should be copied. Remember to create a valid proxy to access data stored in tape.
 2. Use the function Makefile_list(pre_production, basename) to produce a list (not a file!) of filenames contained in the original directory
    stored in the tape. By default, the result is written on a temp local file called filelist_online.txt.
 3. If a selection is necessary in order not to copy all the files of the directory, it must be done before calling the next function.
    The file filelist_online.txt must contain only the filenames to be passed to the next function.
 4. Bring online the files listed in filelist_online.txt calling the function BringOnline(production, tape_filelist="filelist_online.txt") 
 5. FormatFilelistToCopy(production,tape_filelist_tmp, destDir, txtSrmCopyFile) creates a file to be passed to the copy function
 6. SRMCopy(txtSrmCopyFile,decompress_at_destination=False) Is the function to copy the files as written in the list txtSrmCopyFile.
    The optional parameter decompress_at_destination(default value False) allow to unzip copied files if they have the .gz extension.

Finally, files temporally created are deleted and working directory cleaned.

@attention: This function must be executed after initializing a proxy from the UI

@author: N.Tonello
'''

import os
import time
import subprocess
import shutil
import logging
log = logging.getLogger('paudm.tools.grid.storage_utils') 


class ProxyError(Exception):
    pass

class EmptySource(Exception):
    pass


def list_storage_directory(path, config, pattern = '', file_list_file = None):
    
    '''
    Make a list of files to bring online.
       - srmls to list the files available at PIC. This needs a valid proxy. 
         If it is not available, it is automatically initiated, but you should type your personal GRID password.
    @param path: string path to list files
    @param pattern: string pattern to filter i.e. decam or paucam
    @param SRMprefix: string default value "srm://srm.pic.es:8443"
    @param file_list_file: string (optional) filename where to write the srmls result default value
    
    @return: list with dictionaries found with srmls with keywords 'name' and 'path'
    '''
        
    # Use default tmp file or given file (saved)
    if file_list_file == None:
        out_file = "~/tmp_filelist.txt"
    else:
        out_file = file_list_file
    use_srmls = False
    # SRMLS command
    if config['grid']['srm_copy_cmd'] == "srmcp":
        ls_cmd = "srmls "
        use_srmls = True
    elif config['grid']['srm_copy_cmd'] == "lcg-cp":
        ls_cmd = "lcg-ls -D srmv2 -b --connect-timeout=500 "

    if pattern == '':
        cmd = ls_cmd + config['grid']['srm_prefix'] + path + ' | sort -k 2 > ' + out_file
    else:
        cmd = ls_cmd + config['grid']['srm_prefix'] + path + ' | grep ' + pattern + ' | sort -k 2 > ' + out_file
    log.debug("Listing files in %s with %s..." %(path,ls_cmd))
    log.debug(cmd)
    subprocess.call(cmd, shell=True)
    
    lines = os.popen('cat '+ out_file).read().split("\n")
    file_list = []
    for line in lines:

        # Using srmls
        if use_srmls == True:
            try:
                if pattern in line:
                    fits_size = line.split()[0]
                    fits_name = line.split()[1]
                    if fits_size == "" or fits_name == "":
                        continue
                    log.debug("File %s found" %fits_name)
                    if int(fits_size) == 0:
                        log.warning("File %s skipped, it might be corrupted. Filesize is 0."% fits_name)
                        continue
                    file_list.append({'name': os.path.basename(fits_name), 'path':path})  # Example of fits_name:  decam--20-0-g-0.fits.gz
            except Exception:
                log.debug("Could not split this line: %s" %line)
        else:
            try:
                if pattern in line:
                    fits_name = line
                    if fits_name == "":
                        continue
                    log.debug("File %s found" %fits_name)
                    file_list.append({'name': os.path.basename(fits_name), 'path':path})  # Example of fits_name:  decam--20-0-g-0.fits.gz
            except Exception:
                log.debug("Could not split this line: %s" %line)
    if len(file_list) == 0:
        log.warning("No files in storage path")
        
    # Remove previously created out_file, if not given
    if (os.path.exists(out_file)) and (file_list_file == None):
        os.remove(out_file)

    return file_list
    

def bring_online(config_grid, dest_path, file_list, DCAP_prefix = "dcap://dcap.pic.es"):
    '''
    Bring online the files listed in tape_file_list and to be copied.
    This function call the dccp -P command and checks the status of the files every minute.
    It returns only when all listed files are 'ONLINE_AND_NEARLINE'.
    
    @param dest_path: string configuration general IN_PATH
    @param SRM_prefix: string (optional) srm://srm.pic.es:8443
    @param DCAP_prefix: string (optional) gsidcap://gsidcap.pic.es:22128
    @param file_list: dictionary with 'name' as filename and 'path' as storage path with the files to bring online
    '''
    
    log.info("--- PreStage : Bring Online ---")
    
    for file in file_list:
        file_path = file['path']
        # nfs origin
        if file['path'].split('/')[1] == 'nfs':
            shutil.copyfile(os.path.join(file['path'], file['name']), os.path.join(dest_path, file['name']))
            continue
        # pnfs origin
        if file['path'].split('/')[0] == 'pnfs':
            file_path = DCAP_prefix + "/"+file['path']
        elif file['path'].split('/')[1] == 'pnfs':
                file_path = DCAP_prefix + file['path']
        cmd_dcap ='dccp -P ' + os.path.join(file_path, file['name'])
        log.debug(cmd_dcap)
        if subprocess.call(cmd_dcap, shell=True) > 0:
            log.error("System error: Input/output error - \n Check if %s exists."% file['path'])
            return
        
    # Check FilesStatus: they should be all ONLINE at least
    for file in file_list:
        if not os.path.exists(os.path.join(dest_path, file['name'])):
            ready = 'false'
            while ready != 'true':
                if config_grid['srm_copy_cmd'] == "srmcp":
                    srmls_cmd = 'srmls -l ' +  os.path.join(config_grid['srm_prefix'] + file['path'], file['name'])
                    srmls_response = os.popen(srmls_cmd).read()
                    file_status = srmls_response.split()
                    if ('locality:ONLINE_AND_NEARLINE' in file_status) or ('locality:ONLINE' in file_status):
                        ready= 'true'
                        log.info("File %s online (srm)"% file['name'])
                    if ('Error' in file_status):
                        log.error("An error occurred: %s"%srmls_response)
                    else:
                        log.debug("File status: %s"%file_status)
                        time.sleep(30)
                elif config_grid['srm_copy_cmd'] == "lcg-cp":
                    lcg_cmd = 'lcg-ls -l -D srmv2 -b --connect-timeout=500 ' +  os.path.join(config_grid['srm_prefix'] + file['path'], file['name'])
                    lcg_response = os.popen(lcg_cmd).read()
                    file_status = lcg_response.split()
                    if ('ONLINE_AND_NEARLINE' in file_status) or ('ONLINE' in file_status):
                        ready= 'true'
                        log.info("File %s online (lcg)"% file['name'])
                    else:
                        log.debug("File status: %s"%file_status)
                        time.sleep(30)
                else:
                    log.error("An error occurred: %s is not a valid SRM_COPY_CMD"%config_grid['srm_copy_cmd'])
                    break


def srm_download(config_grid, dest_path, file_list = None, decompress_at_destination = False):
    '''
    Copy the files listed with the command specified in common config['grid']['srm_copy_cmd']. 
    
    @param dest_path: NFS destination path
    @param file_list: dictionary with 'name' as filename and 'path' as storage path
    @param decompress_at_destination:  boolean default value False, unzip/unpack the files after copying them
    '''
    
 
    log.debug( "--- Download files (SRM protocol) ---" )
    
    srm_prefix  = config_grid['srm_prefix']
    file_prefix = config_grid['file_prefix']
    if len(file_list) > 0:
        for file in file_list:
            if os.path.exists(os.path.join(dest_path, file['name'])):
                log.info("File %s already in destination path %s"% (file['name'], dest_path))
                continue
            file_tape = os.path.join(srm_prefix + file['path'], file['name'])
            if config_grid['srm_copy_cmd'] == "lcg-cp":
                # pnfs Test
                if dest_path.split('/')[0] == 'pnfs':
                    file_nfs = os.path.join(srm_prefix + dest_path, file['name'])
                elif dest_path.split('/')[1] == 'pnfs':
                    dest_path_list = list(dest_path)[1:]
                    dest_path_formatted = "".join(dest_path_list)
                    file_nfs = os.path.join(srm_prefix + dest_path_formatted, file['name'])
                else:
                    file_nfs  = os.path.join(file_prefix+ dest_path, file['name'])
                srm_cp_cmd = "lcg-cp -b --connect-timeout=500 -D srmv2 "
            elif config_grid['srm_copy_cmd'] == "srmcp":
                # pnfs Test
                if dest_path.split('/')[1] == 'pnfs':
                    file_nfs = os.path.join(srm_prefix + dest_path, file['name'])
                elif dest_path.split('/')[0] == 'pnfs':
                    file_nfs = os.path.join(srm_prefix + "/"+dest_path, file['name'])
                else:
                    file_nfs  = os.path.join(file_prefix+ dest_path, file['name'])                
                srm_cp_cmd = "srmcp"
            if not os.path.exists(os.path.join(dest_path, file['name'])):
                cmd_cp =  srm_cp_cmd + " %s %s" % (file_tape, file_nfs)
                log.debug(cmd_cp)
                if subprocess.call(cmd_cp, shell=True) > 0:
                    error_msg = config_grid['srm_copy_cmd'] + " error while downloading %s to %s" %(file_tape, file_nfs)
                    log.error(error_msg)
                    raise Exception, error_msg
    else:
        log.warning("No files to download. Return")
        return
    
    # Decompress files
    if decompress_at_destination == True :
        log.debug("Uncompressing downloaded fits files..")
        for file in file_list:
            decompress_file(file_name = file['name'], file_path = dest_path)


def decompress_file(file_name, file_path, funpack_bin= "funpack"):
    '''
    Decompress a file located at a given path.
    If gz extension is found, it gunzips the file.
    If fz extension is found, it funpacks the file.
    If the extension is not listed before, ignores it.
    Once completed the decompression, removes the original.
    
    @param file_name: file name
    @param file_path: file path
    @funpack_bin: config value for executable in grid (optional)
    '''
    base, extension = os.path.splitext(file_name)
    
    # Check if already exists the uncompressed version of the file
    if os.path.exists(os.path.join(file_path, base)):
        log.warning("File %s uncompressed at %s. Nothing to do.." %(file_name, file_path))
    else:
        # ZIP file
        if extension == '.gz':
            zipfilename = file_name
            dest_path = file_path
            cmd_unzip = 'gunzip -c ' + os.path.join(dest_path, zipfilename) + ' > ' + dest_path + base
            log.debug(cmd_unzip)
            ret = subprocess.call(cmd_unzip, shell=True)
            if ret == 0:
                log.debug("File %s unzipped. Removing original compressed file...", zipfilename)
                os.remove(os.path.join(dest_path, zipfilename))
                log.debug("Original file removed.")
            else:
                error_msg = "Error while unzipping the file %s with gz" %zipfilename
                log.error(error_msg)
                raise Exception, error_msg
            
        # FPACK file
        elif extension == '.fz':
            # Check if funpack binary code is available
            if not os.path.exists(funpack_bin):
                msg_error = "Could not find funpack code in %s" %funpack_bin
                log.error(msg_error)
                raise Exception, msg_error
                
            cmd_funpack = funpack_bin + ' '+ os.path.join(dest_path, zipfilename)
            log.debug(cmd_funpack) 
            
            if ret == 0:
                log.debug("File %s unpacked. Removing original compressed file...", zipfilename)
                os.remove(os.path.join(dest_path, zipfilename))
                log.debug("Original file removed.")
            else:
                error_msg = "Error while unpacking the file %s with funpack" %zipfilename
                log.error(msg_error)
                raise Exception, msg_error
        else:
            log.debug("Unknown extension for decompression: %s. Ignoring it." %extension)


def download(config_grid, file_list, dest_path = "./", protocol = "SRM"):
    '''
    Function Calls:
        - L{bring_online}
        - L{fts_download}
        - L{srm_download}
    '''
    
    bring_online(config_grid, dest_path, file_list = file_list)
    
    if protocol == 'SRM':
        srm_download(config_grid, dest_path, file_list = file_list, decompress_at_destination=True)
    else:
        msg_error = "Unknown protocol %s" %protocol
        log.error(msg_error)
        raise Exception, msg_error
        
    log.info( 'PreStage DONE' )



# Upload a single file, or all the files of a directory to the storage (srm protocol)
def upload(config, dir_in = None, dir_out= None, file_in = None, file_out = None, transfer_list = None):
    '''
    Upload files from /nfs to storage
    It is possible to specify in the arguments:
        a. file_in and file_out or
        b. file_in and dir_out or
        c. dir_in and dir_out or
        d. transfer_list
        
    parameters:
    @param dir_in : /nfs input directory (where files are)
    @param dir_out : /pnfs output directory (where files will be copied)
    @param file_in : /nfs input file
    @param file_out : /pnfs output file (where it will be copied)
    @param transfer_list : full path of the copyfilelist (format file:////nfs/path/filename srm://...)

    '''
    n_transferred = 0
    n_failed = 0
    list_dir_in = []
    list_failed = []
    
    #When a transfer_list file is given, copy all of the files in the list
    #TBD: list should be passed as a dictionary
    if transfer_list != None:
        if os.path.exists(transfer_list):
            log.info( "Uploading files from list %s ..." % transfer_list)
            str_tmp = "srmcp -copyjobfile="+ transfer_list
            log.debug(str_tmp)
            res = subprocess.call(str_tmp, shell=True)
            if res == 0:
                log.info("...done successfully!")
            else:
                log.error("Something has failed... Sorry!!")
        else:
            n_failed = 1
            log.error("Filelist %s does not exist, cannot upload any file."% transfer_list)

    #When a dir_in is given, copy all the files of the directory
    elif dir_in != None:
        #List files in dir_in, all of them will be transferred
        try:
            list_dir_in = os.listdir(dir_in)
            len_list_in = len(list_dir_in)
            log.info("%d elements found in input directory %s"% (len_list_in, dir_in))

        except Exception:
            log.error( 'Input directory (dir_in) %s not valid. Check it and try again.'% dir_in)

    #When a file_in is given, copy the files to the output directory/file
    elif file_in != None:
        #Check for the existence of the input file
        if os.path.isfile(file_in):
            dir_in = os.path.dirname(file_in)
            element = os.path.basename(file_in)
            list_dir_in = [element]
        else:
            log.error("File_in = %s does not exists.  Check it and try again"% file_in)
            raise Exception

    #Check the validity of the output directory
    if dir_out == None:
        if file_out == None:
            log.error("No destination directory given.\n Please set the file_out or dir_out parameters and try again.")
            list_failed.append(os.path.basename(file_in))
        else:
            dir_out = os.path.dirname(file_out)
            
    if config['grid']['srm_copy_cmd'] == 'lcg-cp':
        cp_cmd = "lcg-cp -b --connect-timeout=500 -D srmv2 "
        ls_cmd = "lcg-ls -b --connect-timeout=500 -D srmv2 --count 10 "
        if dir_out.split("/")[1] == 'pnfs':
            dir_out_list = list(dir_out)[1:]
            dir_out = "".join(dir_out_list)
            
    elif config['grid']['srm_copy_cmd'] == 'srmcp':
        cp_cmd = "srmcp "
        ls_cmd = "srmls "
        if dir_out.split("/")[0] == 'pnfs':
            dir_out = config['grid']['srm_prefix'] + "/" + dir_out
    
    srm_dirout =  config['grid']['srm_prefix']+ dir_out
    if  subprocess.call(ls_cmd + srm_dirout, shell=True) == 1:
        log.warning( "%s does not exist. It will be automatically created."% dir_out)

    #Loop over the files to transfer from /nfs to storage
    for element in list_dir_in:
        if os.path.isdir(os.path.join(dir_in, element)):
            log.info("%s is a directory, do not transfer"% element)
            list_dir_in.remove(element)
        else:
            log.info( "Uploading file %s ..." % element)
            srm_filein = "file:////" + os.path.join(dir_in, element)
            srm_fileout = os.path.join(srm_dirout, element)
            
            str_tmp = cp_cmd + srm_filein + " " + srm_fileout
            log.debug(str_tmp)
            res = subprocess.call(str_tmp, shell=True)
            if res == 0:
                n_transferred = n_transferred + 1
                log.info("...file successfully transferred!")
            elif res ==1:
                log.warning("File %s exists. Check it."% srm_fileout)
                list_failed.append(element)
            else:
                list_failed.append(element)
                log.error( "Upload of file %s has failed... Sorry!!"% element)

                
    log.info("Files successfully transferred: %d / %d "% (n_transferred, len(list_dir_in)))
    log.info("Number of Files not transferred: %d / %d "% (len(list_failed), len(list_dir_in)))

    if n_failed >0:
        raise Exception, "Uploading error for file(s) " + str(list_failed)
    
    log.info("UPLOAD completed!")
    

