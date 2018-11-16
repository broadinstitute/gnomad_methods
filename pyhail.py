#!/usr/bin/env python
from __future__ import print_function
import argparse
import sys
import subprocess
import os
import zipfile
import tempfile

# Great hack for 2.X and 3.X to use input()
try:
    input = raw_input
except NameError:
    pass

try:
    standard_scripts = os.environ['HAIL_SCRIPTS'].split(':')
except Exception:
    standard_scripts = None


def main(args, pass_through_args):
    temp_py = None
    if args.script is None:
        if args.inline is None:
            print('Either --script or --inline is required. Exiting.', file=sys.stderr)
            sys.exit(1)
        if 'print' not in args.inline:
            continue_script = input('No print statement found. Continue? [no] ')
            if not len(continue_script.strip()) or continue_script[0] != 'y':
                sys.exit(1)
        temp_py = tempfile.mkstemp(suffix='.py')
        with open(temp_py[1], 'w') as temp_py_f:
            script = "import hail as hl\nfrom pprint import pprint\n"
            if standard_scripts is not None and any(['gnomad_hail' in x for x in standard_scripts]):
                script = "from gnomad_hail import *\n" + script
            temp_py_f.write(script)
            temp_py_f.write(args.inline)
        script = temp_py[1]
    else:
        script = args.script

    print('Running {} on {}'.format(script, args.cluster))

    hash_string = ''
    try:
        with open(os.path.expanduser(os.environ['HAIL_HASH_LOCATION'])) as f:
            hash_string = f.read().strip()
    except Exception:
        pass

    spark_version = '2.2.0' if args.hail_version == '0.2' else '2.0.2'
    if not hash_string:
        hash_location = 'gs://hail-common/builds/{}/latest-hash/cloudtools-3-spark-{}.txt'.format(args.hail_version, spark_version)
        hash_string = subprocess.check_output(['gsutil', 'cat', hash_location], universal_newlines=True).rstrip()

    if not hash_string:
        print('Could not get hash string', file=sys.stderr)
        sys.exit(1)

    if args.jar is not None:
        jar = args.jar
        jar_file = os.path.basename(jar)
    else:
        jar_file = 'hail-{}-{}-Spark-{}.jar'.format(args.hail_version, hash_string, spark_version)
        jar = 'gs://hail-common/builds/{}/jars/{}'.format(args.hail_version, jar_file)

    all_pyfiles = [args.zip if args.zip is not None else 'gs://hail-common/builds/{0}/python/hail-{0}-{1}.zip'.format(args.hail_version, hash_string)]

    pyfiles = []
    if args.add_scripts:
        pyfiles.extend([os.path.expanduser(x) for x in args.add_scripts.split(',')])
    if standard_scripts is not None:
        pyfiles.extend(standard_scripts)
    if pyfiles:
        tfile = tempfile.mkstemp(suffix='.zip', prefix='pyscripts_')[1]
        print(tfile)
        zipf = zipfile.ZipFile(tfile, 'w', zipfile.ZIP_DEFLATED)
        for hail_script_entry in pyfiles:
            if hail_script_entry.endswith('.py'):
                zipf.write(hail_script_entry, arcname=os.path.basename(hail_script_entry))
            else:
                for root, _, files in os.walk(hail_script_entry):
                    for pyfile in files:
                        if pyfile.endswith('.py'):
                            zipf.write(os.path.join(root, pyfile),
                                       os.path.relpath(os.path.join(root, pyfile),
                                                       os.path.join(hail_script_entry, '..')))
        zipf.close()
        all_pyfiles.append(tfile)

    print('Using JAR: {} and files:\n{}'.format(jar, '\n'.join(pyfiles)))

    job = ['gcloud', 'dataproc', 'jobs', 'submit', 'pyspark', script,
           '--cluster', args.cluster,
           '--files={}'.format(jar),
           '--py-files={}'.format(','.join(all_pyfiles)),
           # '--properties={}'.format(','.join(spark_properties)),
           '--driver-log-levels', 'root=FATAL,is.hail=INFO'
           ]
    spark_properties = []
    spark_properties.extend(['spark.{}=./{}'.format(x, jar_file) for x in ('executor.extraClassPath', 'driver.extraClassPath', 'files')])
    spark_properties.append('spark.submit.pyFiles=./{}'.format(all_pyfiles[0]))
    if args.spark_conf:
        spark_properties.extend(args.spark_conf.split(','))
    job.append('--properties={}'.format(','.join(spark_properties)))

    if pass_through_args is not None:
        job.append('--')
        job.extend(pass_through_args)
    print('Running: {}'.format(' '.join(job)))

    subprocess.check_output(job)
    if temp_py is not None:
        os.remove(temp_py[1])


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--script', '--input', '-i', help='Script to run')
    parser.add_argument('--inline', help='Inline script to run')
    parser.add_argument('--cluster', help='Which cluster to run on', required=True)

    hail_script_options = parser.add_argument_group('Additional hail script options')
    hail_script_options.add_argument('--hail_version', help='Hail version (0.1 or 0.2)', default='0.2')
    hail_script_options.add_argument('--jar', help='Jar file to use')
    hail_script_options.add_argument('--zip', help='Hail zip file to use')
    hail_script_options.add_argument('--add_scripts', help='Comma-separated list of additional python scripts to add.')
    hail_script_options.add_argument('--spark_conf', help='Comma-separated list of additional spark configurations to pass.')
    args, pass_through_args = parser.parse_known_args()
    main(args, pass_through_args)
