#####################################################################################################################
##
## The module is part of SiteProf2 framework, for processing structures& sequences of biological macromolecules.
## It may not be copied without explicit permission of the author.
## Copyright: Grzegorz M. Koczyk (2005-2007)
##
#####################################################################################################################
from __future__ import with_statement
import sys, os, os.path
############################################################################################
# Dummy object corresponding to C struct construct
class StructObject(object):
    def __init__(self, **dictargs):    
        for k, v in dictargs.items():
            object.__setattr__(self,k, v)
            
    def __setattr__(self, k, v):
        object.__setattr__(self, k, v)
    def __str__(self):
        return "( " + ", ".join(["%s: %s" % (k, self.__dict__[k]) for k in sorted(self.__dict__)]) + " )";
    def __repr__(self):
        return self.__str__()
    
def handlify(mode='r', close=False, is_method=False):
    ''' Decorator on a parse/read function.method, so that if a string is supplied it will attempt to open a file
    and use it instead '''
    # Intermediate function is needed so that decorator creation can be customised, by options above
    def intermediate(func):
        if is_method:
            def result(self, fh, *args, **kwargs):
                if isinstance(fh, str):
                    with open(fh, mode, newline='') as fh:
                        r = func(self, fh, *args, **kwargs)
                else:
                    r = func(self, fh, *args, **kwargs)
                return r
            return result
        else:                
            def result(fh, *args, **kwargs):
                if isinstance(fh, str):
                    with open(fh, mode, newline='') as fh:
                        r = func(fh, *args, **kwargs)
                else:
                    r = func(fh, *args, **kwargs)
                return r
        return result
    return intermediate

from collections import defaultdict
def innerDefaultdict(t=None):
    def __inner():
        return defaultdict(t)
    return __inner

from contextlib import contextmanager
@contextmanager
def changedDir(new_dir):
    try:
        old_dir = os.getcwd()
        os.chdir(new_dir)
        yield new_dir
    finally:
        os.chdir(old_dir)

import tempfile, shutil
@contextmanager
def tempDir(suffix='.dhcl', prefix='tmp.', dir='.', change_dir=True):
    try:
        old_dir = os.getcwd()
        tmp_dir = tempfile.mkdtemp(suffix, prefix, dir)
        if change_dir:
            os.chdir(tmp_dir)
        yield tmp_dir
    finally:
        if change_dir:
            os.chdir(old_dir)
            shutil.rmtree(tmp_dir)

# Helper function to provide a temporary directory to a function
import tempfile, shutil, types
import contextlib
@contextlib.contextmanager
def tmpfile(give_name=False, ext='.tmp', prefix='tmp', dir=None):
    try:
        tmp_fh, tmp_fname = tempfile.mkstemp(suffix=ext,dir=dir, prefix=prefix)
        if give_name:
            yield tmp_fh, tmp_fname
        else:
            yield tmp_fh
    finally:
        try:
            os.close(tmp_fh)
        except:
            pass
        try:
            os.unlink(tmp_fname)
        except:
            pass

@contextlib.contextmanager
def tmpdir(suffix='.tmp', prefix='tmp', dir='/tmp'):
    try:
        tmp_dir = tempfile.mkdtemp(suffix=suffix, prefix=prefix, dir=dir)
        yield tmp_dir
    finally:
        try:
            shutil.rmtree(tmp_dir)
        except:
            pass

def provideTempdir(suffix='.tmp', prefix='tmp', dir='/tmp', delete_dir=True, change_dir=False, is_method=False, is_generator=False):
    def intermediate(func):
        if not is_method:
            if not is_generator:
                def result(*args, **kwargs):
                    dirname = tempfile.mkdtemp(suffix=suffix, prefix=prefix, dir=dir)
                    old_dir = os.getcwd()
                    try:
                        if change_dir:
                            os.chdir(dirname)
                        return func(dirname, *args, **kwargs)
                    finally:
                        if change_dir:
                            os.chdir(old_dir)
                        if delete_dir and os.path.exists(dirname):
                            shutil.rmtree(dirname)
            else:
                def result(*args, **kwargs):
                    dirname = tempfile.mkdtemp(suffix=suffix, prefix=prefix, dir=dir)
                    old_dir = os.getcwd()
                    try:
                        if change_dir:
                            os.chdir(dirname)
                        for r in func(dirname, *args, **kwargs):
                            yield r
                        if change_dir:
                            os.chdir(old_dir)
                    finally:
                        if change_dir:
                            os.chdir(old_dir)
                        if delete_dir and os.path.exists(dirname):
                            shutil.rmtree(dirname)
        else:
            if not is_generator:
                def result(self, *args, **kwargs):
                    dirname = tempfile.mkdtemp(suffix=suffix, prefix=prefix, dir=dir)
                    old_dir = os.getcwd()
                    try:
                        if change_dir:
                            os.chdir(dirname)
                        return func(self, dirname, *args, **kwargs)
                    finally:
                        if change_dir:
                            os.chdir(old_dir)
                        if delete_dir and os.path.exists(dirname):
                            shutil.rmtree(dirname)
            else:
                def result(self, *args, **kwargs):
                    dirname = tempfile.mkdtemp(suffix=suffix, prefix=prefix, dir=dir)
                    old_dir = os.getcwd()
                    try:
                        if change_dir:
                            os.chdir(dirname)
                        for r in func(self, dirname, *args, **kwargs):
                            yield r
                        if change_dir:
                            os.chdir(old_dir)
                    finally:
                        if change_dir:
                            os.chdir(old_dir)
                        if delete_dir and os.path.exists(dirname):
                            shutil.rmtree(dirname)
        return result
    return intermediate

def provideTempfile(give_name=False, ext='.tmp', prefix='tmp', dir=None, is_method=False):
    def intermediate(func):
        if not is_method:
            def result(*args, **kwargs):
                fh, fname = tempfile.mkstemp(suffix=ext,dir=dir, prefix=prefix)
                try:
                    if give_name:
                        return func(fh, fname, *args, **kwargs)
                    else:
                        return func(fh, *args, **kwargs)
                finally:
                    try:
                        os.close(fh)
                    except:            
                        pass
                    try:
                        os.unlink(fname)
                    except:
                        pass
        else:
            def result(self, *args, **kwargs):
                fh, fname = tempfile.mkstemp(suffix=ext,dir=dir, prefix=prefix)
                try:
                    if give_name:
                        return func(self, fh, fname, *args, **kwargs)
                    else:
                        return func(self, fh, *args, **kwargs)
                finally:
                    try:
                        os.close(fh)
                    except:            
                        pass
                    try:
                        os.unlink(fname)
                    except:
                        pass
        return result
    return intermediate
         
#######################################################################################################################
## Parsing/writing basic forms of list, set and dictionary (1/line, dicts in tab/comma/whatever-separated files)
#######################################################################################################################
@handlify(is_method=False)
def parseSimpleList(fh, type_=str):
    result = list()
    for line in fh:
        line = line.strip()
        if line:
            result.append( type_(line) )
    return result

@handlify(mode='w', is_method=False)
def writeSimpleList(fh, lyst):
    for v in lyst:
        print >> fh, "%s" % (v)

@handlify(is_method=False)
def parseSimpleSet(fh):
    result = set()
    for line in fh:
        line = line.strip()
        if line:
            result.add(line)
    return result

@handlify(is_method=False)
def parseSimpleDict(fh, separator='\t', key_f=None, value_f=None):
    ''' Parses a dictionary in form of multiple character-separated values (one line, one key-value pair)'''    
    result = {}
    for line in fh:
        line = line.strip()
        fs = line.split(separator, 1)
        k, v = fs[0], fs[1]
        if key_f is not None:
            k = key_f(k)
        if value_f is not None:
            v = value_f(v)
        result[k] = v
    return result

@handlify(mode='w', is_method=False)
def writeSimpleDict(fh, dikt, separator='\t', **kwds):
    ''' Writes a dictionary in form of multiple character-separated values (one line, one key-value pair)'''
    for k,v in dikt.iteritems():
        print >> fh, "%s%s%s" % (k,separator,v)
    if kwds.has_key('eof'):
        print >> fh, kwds['eof'] 


def iterateDir(in_dir, prefix=None, suffix=None):
    ''' Process a directory into a set of filenames (optionally by suffix or prefix) '''
    for fname in os.listdir(in_dir):
        if (prefix is None or fname.startswith(prefix)) and (suffix is None or fname.endswith(suffix)):
            full_fname = os.path.abspath(os.path.join(in_dir, fname))
            yield full_fname

@contextmanager
def condMakeDir(out_dir):
    out_dir = os.path.abspath(out_dir)
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    with changedDir(out_dir):
        yield out_dir
        
def fname2pdb(fname):
    return os.path.basename(fname).split('.')[0].lower()

def stringize(it):
    ''' Convert an iterables\' elements to strings. '''
    for i in it:
        yield str(i)

@handlify(mode='r', is_method=False)
def endsWithDone(fh):
    done = False
    for line in fh:
        if line.endswith('#DONE'):
            done=True
            break
    return done

def existsDone(fname):
    ''' Does a .done marker file exists '''
    return os.path.exists(fname+'.done')

def markAsDone(fname):
    '''  Creates a .done marker file '''
    open(fname+'.done', 'w').close()

def hasContent(fname):
    return (os.path.exists(fname) and os.path.isfile(fname) and not os.path.getsize(fname)==0)

def parseResidueId(s):
    ''' Parse residue identifier into (number, insertion code) tuple.'''
    s = s.strip()
    if s[-1] in '0123456789':
        resno, icode = int(s), ' '
    else:
        resno, icode = int(s[:-1]), s[-1]
    return resno,icode

#######################################################################################################################
## Parsing/writing basic forms of list, set and dictionary (1/line, dicts in tab/comma/whatever-separated files)
#######################################################################################################################
@handlify(is_method=False)
def parseSimpleList(fh, type_=str):
    result = list()
    for line in fh:
        line = line.strip()
        if line:
            result.append( type_(line) )
    return result

@handlify(mode='w', is_method=False)
def writeSimpleList(fh, lyst):
    for v in lyst:
        print >> fh, "%s" % (v)

@handlify(is_method=False)
def parseSimpleSet(fh):
    result = set()
    for line in fh:
        line = line.strip()
        if line:
            result.add(line)
    return result

@handlify(is_method=False)
def parseSimpleDict(fh, separator='\t', key_f=None, value_f=None):
    ''' Parses a dictionary in form of multiple character-separated values (one line, one key-value pair)'''    
    result = {}
    for line in fh:
        line = line.strip()
        fs = line.split(separator, 1)
        k, v = fs[0], fs[1]
        if key_f is not None:
            k = key_f(k)
        if value_f is not None:
            v = value_f(v)
        result[k] = v
    return result

@handlify(mode='w', is_method=False)
def writeSimpleDict(fh, dikt, separator='\t', **kwds):
    ''' Writes a dictionary in form of multiple character-separated values (one line, one key-value pair)'''
    for k,v in dikt.iteritems():
        print >> fh, "%s%s%s" % (k,separator,v)
    if kwds.has_key('eof'):
        print >> fh, kwds['eof'] 

#######################################################################################################################
## IniOptionParser - extended ConfigParser defaulting to .ini file for values
#######################################################################################################################
import configparser, optparse

class IniOptionParser(optparse.OptionParser, object):
    ''' Subclassed OptionParser from optparse, which can take arguments from the parsed .INI config file. '''
    def __init__(self, inifile, *args, **kwargs):
        optparse.OptionParser.__init__(self, *args, **kwargs)
        self.add_option("--inifile", dest="inifile",
                        default=inifile,
                        help="specify a file containing further specified options (if invalid print a warning and proceed) [%s]" % (inifile)
                        )
    def _dictify(self, ini_opts, types={}):
        result = {}
        for section in ini_opts.sections():
            if section=='__MAIN__':
                inner = result
            else:
                inner = result[section] = {}
            for option, value in ini_opts.items(section):
                if section in types and option in types[section]:
                    inner[option] = types[section][option](value)
                else:
                    inner[option] = value
        return result
    def parse_args(self, types={}):
        (options, args) = super(IniOptionParser, self).parse_args()
        if os.path.exists(options.inifile):
            parser = configparser.ConfigParser()
            parser.read([options.inifile])
            ini_options = self._dictify(parser, types)
            for key, value in ini_options.iteritems():
                self.defaults[key] = value
            (options, args) = super(IniOptionParser, self).parse_args()
        return (options, args)

def makeListTokenizer(separator=",", create_f=str):
    def inner(data):
        return [create_f(x) for x in data.split(separator)]
    return inner

def separatingCallback(option, opt_str, value, parser, create_f=str, separator=","):
    if value:
        setattr(parser.values, option.dest, makeListTokenizer(separator, create_f)(value))

EXTENSIONS = {
    'vdw': '.residue.cs',
    'loops': '.loops',
    'locks': '.locks',
    'domains': '.domains'
    }
 
