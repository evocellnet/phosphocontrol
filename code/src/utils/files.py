"""
Functions to deal with the filesystem
"""
import os

def is_file_empty(file_path):
    return (os.stat(str(file_path)).st_size == 0)
