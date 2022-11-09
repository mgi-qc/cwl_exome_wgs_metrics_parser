import os, csv, subprocess, datetime, argparse, glob, sys

# argument input
desc_str = """
        Program to parse cwl metrics.
    """
parser = argparse.ArgumentParser(description=desc_str)
parser.add_argument("-w", type=str, help='woid')
parser.add_argument("-fw", type=str, help='file of woid\'s (without header)')
parser.add_argument("-a", type=str, help='anp')
parser.add_argument("-fa", type=str, help='file of AnP\'s (without header)')
parser.add_argument("-m", type=str, help='model group id')
parser.add_argument("-fm", type=str, help='file of model group id\'s')
parser.add_argument("-l", help='create library output file', action='store_true')
parser.add_argument("-e", help='Run exome report maker', action='store_true')
parser.add_argument("-wgs", help='Run wgs report maker', action='store_true')
args = parser.parse_args()

id_list = []

# assign woid
if not args.w and not args.a and not args.fw and not args.fa and not args.m and not args.fm and not args.l:
    sys.exit('usage:\npython3.5 cwlbleau.py -w <woid>\npython3.5 cwlbleau.py -f <file of woids>')

elif args.w:
    anp_or_woid = "WorkOrder"
    id_list.append(args.w)
elif args.a:
    anp_or_woid = "AnP"
    id_list.append(args.a)
elif args.fw:
    anp_or_woid = "WorkOrder"
    with open(args.fw, 'r') as infilecsv:
        for line in infilecsv:
            id_list.append(line.rstrip())
elif args.fa:
    anp_or_woid = "AnP"
    with open(args.fa) as infilecsv:
        for line in infilecsv:
            id_list.append(line.rstrip())
elif args.m:
    anp_or_woid = "Model_Group_ID"
    id_list.append(args.m)
elif args.fm:
    anp_or_woid = "Model_Group_ID"
    with open(args.fm, 'r') as infilecsv:
        for line in infilecsv:
            id_list.append(line.rstrip())

# set working dir, results dic, date
working_dir = os.getcwd()
mm_dd_yy = datetime.datetime.now().strftime("%m%d%y")
mmddyy_slash = datetime.datetime.now().strftime("%m/%d/%y")


# Fucntions pull metrics from last succeeded build dir
def verify_bamid(infile):
    with open(infile, 'r') as infilecsv:
        infile_reader = csv.DictReader(infilecsv, delimiter='\t')
        for line in infile_reader:
            results['SEQ_ID'] = line['#SEQ_ID']
            results['FREEMIX'] = line['FREEMIX']
    return results


def insert_size_metrics(infile):
    count = 0
    with open(infile, 'r') as infilecsv:
        infile_reader = csv.reader(infilecsv, delimiter='\t')
        for line in infile_reader:
            if 'MEAN_INSERT_SIZE' in line:
                header_position = line.index("MEAN_INSERT_SIZE")
            if 'MEAN_INSERT_SIZE' in line and count == 0:
                data = (next(infile_reader))
                results['MEAN_INSERT_SIZE'] = data[header_position]
                results['STANDARD_DEVIATION'] = data[header_position+1]
                count = 1
    return results


def flagstat_out(infile):
    with open(infile, 'r') as infilecsv:

        infile_reader = csv.reader(infilecsv)

        mapped_rate = float()
        mapped_int = float()
        properly_paired_rate = float()
        wmmtadc = float()

        for line in infile_reader:

            if 'mapped' in line[0] and '% : N/A' in line[0]:
                mapped_line = line[0]
                mapped_split = mapped_line.split('(')
                mapped_rate = float(mapped_split[1].split(':')[0].strip()[:-1])
                mapped_int = float(mapped_split[0].split('+')[0].strip())
                results['mapped_rate'] = mapped_rate

            if 'properly paired' in line[0]:
                properly_paired_line = line[0]
                properly_paired_split = properly_paired_line.split('(')
                properly_paired_rate = float(properly_paired_split[1].split(':')[0].strip()[:-1])
                results['properly_paired-rate'] = properly_paired_rate

            if 'with mate mapped to a different chr' in line[0] and 'mapQ>=' not in line[0]:
                wmmtadc_line = line[0]
                wmmtadc = float(wmmtadc_line.split('+')[0].strip())

        discordant_rate = mapped_rate - properly_paired_rate
        results['discordant_rate'] = discordant_rate

        inter_chromosomal_pairing_rate = wmmtadc / mapped_int
        results['inter-chromosomal_Pairing rate'] = inter_chromosomal_pairing_rate

    return results


def mark_dups_metrics(infile, sample):
    count = 0
    percent_total = 0
    with open(infile, 'r') as infilecsv:
        infile_reader = csv.reader(infilecsv, delimiter='\t')
        percent_duplication = float()
        for line in infile_reader:
            if line and sample in line[0]:
                percent_total += float(line[8])
                count += 1

        results['PERCENT_DUPLICATION'] = (percent_total/count)
        percent_duplication = (percent_total/count)

    return percent_duplication


def gcbias_metrics_summary(infile):
    with open(infile, 'r') as infilecsv:
        infile_reader = csv.reader(infilecsv, delimiter='\t')
        for line in infile_reader:
            if 'ALIGNED_READS' in line:
                data = next(infile_reader)
                results['ALIGNED_READS'] = data[4]
    return results


def alignment_summary_metrics(infile):
    first_count = 0
    second_count = 0
    pair_count = 0
    with open(infile, 'r') as infilecsv:
        infile_reader = csv.reader(infilecsv, delimiter='\t')
        pf_aligned_bases = float()
        for line in infile_reader:
            if 'FIRST_OF_PAIR' in line and first_count == 0:
                results['FOP: PF_MISMATCH_RATE'] = line[12]
                first_count = 1
            if 'SECOND_OF_PAIR' in line and second_count == 0:
                results['SOP: PF_MISMATCH_RATE'] = line[12]
                second_count = 1
            if 'PAIR' in line and not '_' in line and pair_count == 0:
                pf_aligned_bases = int(line[7])
                results['TOTAL_READS'] = line[1]
                results['PF_READS'] = line[2]
                results['PF_READS_ALIGNED'] = line[5]
                results['PF_ALIGNED_BASES'] = line[7]
                results['PF_HQ_ALIGNED_Q20_BASE'] = line[10]
                results['PCT_ADAPTER'] = line[23]
                pair_count = 1
    return pf_aligned_bases


def wgs_metrics(infile):
    with open(infile, 'r') as infilecsv:
        infile_reader = csv.reader(infilecsv, delimiter='\t')
        genome_territory = int()
        for line in infile_reader:
            if 'GENOME_TERRITORY' in line:
                data = next(infile_reader)
                genome_territory = int(data[0])
                results['GENOME_TERRITORY'] = data[0]
                results['MEAN_COVERAGE'] = data[1]
                results['SD_COVERAGE'] = data[2]
                results['PCT_10X'] = data[14]
                results['PCT_20X'] = data[16]
                results['PCT_30X'] = data[18]
                results['HET_SNP_SENSITIVITY'] = data[26]
                results['HET_SNP_Q'] = data[27]
    return genome_territory


def hs_metrics(infile):
    with open(infile, 'r') as infilecsv:
        infile_reader = csv.reader(infilecsv, delimiter='\t')
        for line in infile_reader:
            if 'BAIT_SET' in line:
                hs_metrics_header = line
                # hs_metrics_data_one = next(infile_reader)
                hs_metrics_data_two = next(infile_reader)
                hs_metrics_dict = dict(zip(hs_metrics_header, hs_metrics_data_two))
        for metric in hs_metrics_dict:
            results[metric] = hs_metrics_dict[metric]
    return hs_metrics_header


def write_results(results_dict, outfile, header_list):
    if not os.path.isfile(outfile):
        with open(cwd_metrics_outfile, 'w') as outfilecsv:
            metrics_writer = csv.DictWriter(outfilecsv, fieldnames=header_list, delimiter='\t')
            metrics_writer.writeheader()
            metrics_writer.writerow(results_dict)
    elif os.path.isfile(outfile):
        with open(cwd_metrics_outfile, 'a') as outfilecsv:
            metrics_writer = csv.DictWriter(outfilecsv, fieldnames=header_list, delimiter='\t')
            metrics_writer.writerow(results_dict)
    return


def write_library_results(results_dict, outfile, header_list):
    if 'Library' not in header_list[0]:
        header_list.insert(0, 'Library')
    if not os.path.isfile(outfile):
        with open(outfile, 'w') as outfilecsv:
            metrics_writer = csv.DictWriter(outfilecsv, fieldnames=header_list, delimiter='\t')
            metrics_writer.writeheader()
            metrics_writer.writerow(results_dict)
    elif os.path.isfile(outfile):
        with open(outfile, 'a') as outfilecsv:
            metrics_writer = csv.DictWriter(outfilecsv, fieldnames=header_list, delimiter='\t')
            metrics_writer.writerow(results_dict)
    return


def file_pick(list_of_files, file_type):

    list_of_files = [file.split('/')[-1] for file in list_of_files]
    print('\n----------\n')
    print('There are {} {} files, please select a file:'.format(len(list_of_files), file_type))
    for number, file in enumerate(list_of_files, 1):
        print('{}. {}'.format(number, file))
    print('\nPlease select a file:')

    while True:
        # print files to terminal, with no, have user pick file, return file name
        user_input = input()
        if user_input.isdigit() and int(user_input) > 0 and int(user_input) - 1 < len(list_of_files):
            filename_list = list_of_files[int(user_input) - 1].split('.')
            file_flag = [i for i, s in enumerate(filename_list) if file_type in s]
            return filename_list[file_flag[0]]

        print('Please pick a file number between 1 and {}:'.format(len(list_of_files)))


def file_check(file_dir):
    exome_qc_files_dict = {}
    exome_qc_file_list = ['HsMetrics', 'mark_dups_metrics', 'GcBiasMetricsSummary',
                          'AlignmentSummaryMetrics', 'WgsMetrics', 'InsertSizeMetrics',
                          'VerifyBamId.selfSM', 'flagstat', 'cram']

    while True:
        for qc_file in exome_qc_file_list:
            found_files = glob.glob(file_dir+'*{}*'.format(qc_file))

        # check each file used,
            if not found_files:
                print('\n----------\n')
                print('{} file not found, result will be FNF.'.format(qc_file))
                exome_qc_files_dict[qc_file] = qc_file
                continue

            if len(found_files) > 1:
                exome_qc_files_dict[qc_file] = file_pick(found_files, qc_file)
            else:
                print('\n----------\n')
                print('Only 1 {file_name} file found. QC\'ing with {file_name}.'.format(file_name=qc_file))
                filename = found_files[0].split('/')[-1]

                if qc_file == 'VerifyBamId.selfSM':
                    exome_qc_files_dict[qc_file] = 'VerifyBamId.selfSM'
                    continue

                filename_list = filename.split('.')
                file_flag = [i for i, s in enumerate(filename_list) if qc_file in s]
                exome_qc_files_dict[qc_file] = filename_list[file_flag[0]]

        # if there's only one, assign to file dict
        # else, check with user which file type to use
        # return file profile to qc, user accepts.
        print('\n**********')
        print('QC Profile\n')
        for file_name, file in sorted(exome_qc_files_dict.items()):
            print('{}: {}'.format(file_name, file))
        file_confirm = input('\nconfirm file selection y, anything else to re-start:\n')

        if 'y' in file_confirm.lower():
            print('**********')
            return exome_qc_files_dict


def get_modelgroup_results_dir(model_list):
    dir_found = False
    for model in model_list:
        results_dir = model.split('\t')[3]
        if os.path.isdir(os.path.join(model.split('\t')[3], 'results/')):
            dir_found = True
            return os.path.join(model.split('\t')[3], 'results/')
    if not dir_found:
        sys.exit('No Results directories found for any model groups.')


# Aye ordered these columns
met_wgs_header = ['Admin', 'WorkOrder', 'date_QC', 'sample_name', 'common_name', 'model_name', 'last_succeeded_build',
                  'data_directory',
                  'cram_file', 'status', 'ALIGNED_READS', 'mapped_rate', 'FOP: PF_MISMATCH_RATE',
                  'SOP: PF_MISMATCH_RATE',
                  'FREEMIX', 'HAPLOID COVERAGE', 'PCT_10X', 'PCT_20X', 'PCT_30X', 'discordant_rate',
                  'inter-chromosomal_Pairing rate', 'HET_SNP_Q', 'HET_SNP_SENSITIVITY',
                  'MEAN_COVERAGE', 'SD_COVERAGE', 'MEAN_INSERT_SIZE', 'STANDARD_DEVIATION', 'PCT_ADAPTER', 'PF_READS',
                  'PF_ALIGNED_BASES', 'PERCENT_DUPLICATION', 'TOTAL_READS', 'properly_paired-rate',
                  'PF_HQ_ALIGNED_Q20_BASE', 'PF_READS_ALIGNED', 'GENOME_TERRITORY', 'SEQ_ID']

met_wgs_header[1] = anp_or_woid


for id in id_list:

    if args.w or args.fw:
        print('\ncwlbleau\'ing: {}'.format(id))
        model_groups_id = 'model_groups.project.id=' + id

        # run genome model command to generate model info
        model_groups = subprocess.check_output(['genome', 'model', 'list', model_groups_id, "--show",
                                                "last_succeeded_build.id,name,status,last_succeeded_build.data_directory,"
                                                "subject.name,subject.common_name", "--style=tsv", "--nohead"]).decode('utf-8').splitlines()

    if args.a or args.fa:
        print('\ncwlbleau\'ing: {}'.format(id))
        model_groups_id = 'analysis_project.id=' + id

        # run genome model command to generate model info
        model_groups = subprocess.check_output(['genome', 'model', 'list', model_groups_id, "--show",
                                                "last_succeeded_build.id,name,status,last_succeeded_build.data_directory,"
                                                "subject.name,subject.common_name ", "--style=tsv", "--nohead"]).decode('utf-8').splitlines()

    if args.m or args.fm:
        print('\ncwlbleau\'ing: {}'.format(id))
        model_groups_id = 'model_groups.id=' + id

        # run genome model command to generate model info
        model_groups = subprocess.check_output(['genome', 'model', 'list', model_groups_id, "--show",
                                                "last_succeeded_build.id,name,status,last_succeeded_build.data_directory,"
                                                "subject.name,subject.common_name ", "--style=tsv", "--nohead"]).decode('utf-8').splitlines()

    model_group_dir = get_modelgroup_results_dir(model_groups)

    if not os.path.isdir(model_group_dir):
        print('SKIPPING EXOME QC FOR: {}'.format(id))
        print('{}/results/ not found.'.format(model_group_dir))
        print('\n----------')
        continue
    else:
        print('Using {} directory for exome file generation.'.format(model_group_dir))
    exome_qc_files = file_check(model_group_dir)

    # create outfiles
    cwd_metrics_outfile = id + '.cwl.metrics.' + mm_dd_yy + '.tsv'
    print('\nMetrics outfile: {}\n'.format(cwd_metrics_outfile))
    print('----------')
    if os.path.isfile(cwd_metrics_outfile):
        os.remove(cwd_metrics_outfile)

    cwd_metrics_library_outfile = id + '.cwl.metrics.library.' + mm_dd_yy + '.tsv'
    if os.path.isfile(cwd_metrics_library_outfile):
        os.remove(cwd_metrics_library_outfile)

    # Admin name
    # disapled admin query because of lims/gms docker container issue
    ap_new = "NA"
    if args.w or args.fw:
        admin_collections = subprocess.check_output(["wo_info", "--report", "billing", "--woid", id]).decode(
            'utf-8').splitlines()
        for ap in admin_collections:
            if 'Administration Project' in ap:
                ap_new = ap.split(':')[1].strip()

    # call methods to generate results
    for line in model_groups:

        results = {}

        info = line.split('\t')

        if 'Succeeded' in info[2]:

            results.clear()
            results['Admin'] = ap_new
            results[anp_or_woid] = id
            results['date_QC'] = mmddyy_slash
            results['last_succeeded_build'] = info[0]
            results['model_name'] = info[1]
            results['status'] = info[2]
            results['data_directory'] = info[3]
            results['sample_name'] = info[4]
            results['common_name'] = info[5]
            if '<NULL>' in info[5]:
                results['common_name'] = 'NA'

            if not os.path.isdir(info[3] + '/results'):
                print('{} directory not found'.format(info[3] + '/results'))
                continue

            os.chdir(info[3] + '/results')

            results['cram_file'] = 'NA'
            cram_file = glob.glob('*{}'.format(exome_qc_files['cram']))
            if not cram_file:
                cram_file = 'NA'
            else:
                cram_file = cram_file[0]
            if os.path.isfile(cram_file):
                results['cram_file'] = os.getcwd() + '/{}'.format(cram_file)

            vbi = glob.glob('*{}*'.format(exome_qc_files['VerifyBamId.selfSM']))
            if not vbi:
                vbi = exome_qc_files['VerifyBamId.selfSM']
            else:
                vbi = vbi[0]
            if os.path.isfile(vbi):
                verify_bamid(vbi)
            else:
                results['SEQ_ID'] = 'FNF'
                results['FREEMIX'] = 'FNF'

            ism = glob.glob('*{}*'.format(exome_qc_files['InsertSizeMetrics']))
            if not ism:
                ism = exome_qc_files['InsertSizeMetrics']
            else:
                ism = ism[0]
            if os.path.isfile(exome_qc_files['InsertSizeMetrics']):
                insert_size_metrics(exome_qc_files['InsertSizeMetrics'])
            else:
                results['MEAN_INSERT_SIZE'] = 'FNF'
                results['STANDARD_DEVIATION'] = 'FNF'

            flagstat = glob.glob('*{}*'.format(exome_qc_files['flagstat']))
            if not flagstat:
                flagstat = exome_qc_files['flagstat']
            else:
                flagstat = flagstat[0]
            if os.path.isfile(flagstat):
                flagstat_out(flagstat)
            else:
                results['mapped_rate'] = 'FNF'
                results['properly_paired-rate'] = 'FNF'
                results['discordant_rate'] = 'FNF'
                results['inter-chromosomal_Pairing rate'] = 'FNF'

            mdp = glob.glob('*{}*'.format(exome_qc_files['mark_dups_metrics']))
            if not mdp:
                mdp = exome_qc_files['mark_dups_metrics']
            else:
                mdp = mdp[0]
            if os.path.isfile(mdp):
                perc_dup = mark_dups_metrics(mdp, info[4])
            else:
                results['PERCENT_DUPLICATION'] = 'FNF'
                perc_dup = False

            gbms = glob.glob('*{}*'.format(exome_qc_files['GcBiasMetricsSummary']))
            if not gbms:
                gbms = exome_qc_files['GcBiasMetricsSummary']
            else:
                gbms = gbms[0]
            if os.path.isfile(gbms):
                gcbias_metrics_summary(gbms)
            else:
                results['ALIGNED_READS'] = 'FNF'

            asm = glob.glob('*{}*'.format(exome_qc_files['AlignmentSummaryMetrics']))
            if not asm:
                asm = exome_qc_files['AlignmentSummaryMetrics']
            else:
                asm = asm[0]
            if os.path.isfile(asm):
                pfalgnbases = alignment_summary_metrics(asm)
            else:
                results['TOTAL_READS'] = 'FNF'
                results['PF_READS'] = 'FNF'
                results['PF_READS_ALIGNED'] = 'FNF'
                results['PF_ALIGNED_BASES'] = 'FNF'
                results['PF_HQ_ALIGNED_Q20_BASE'] = 'FNF'
                results['PCT_ADAPTER'] = 'FNF'
                pfalgnbases = False

            wm = glob.glob('*{}*'.format(exome_qc_files['WgsMetrics']))
            if not wm:
                wm = exome_qc_files['WgsMetrics']
            else:
                wm = wm[0]
            if os.path.isfile(wm):
                genome_terr = wgs_metrics(wm)
            else:
                results['GENOME_TERRITORY'] = 'FNF'
                results['MEAN_COVERAGE'] = 'FNF'
                results['SD_COVERAGE'] = 'FNF'
                results['PCT_10X'] = 'FNF'
                results['PCT_20X'] = 'FNF'
                results['PCT_30X'] = 'FNF'
                results['HET_SNP_SENSITIVITY'] = 'FNF'
                results['HET_SNP_Q'] = 'FNF'
                genome_terr = False

            if pfalgnbases and perc_dup and genome_terr:
                haploid_coverage = pfalgnbases * ((1 - perc_dup)/genome_terr)
                results['HAPLOID COVERAGE'] = haploid_coverage
            else:
                results['HAPLOID COVERAGE'] = 'FNF'

            hsm = glob.glob('*{}*'.format(exome_qc_files['HsMetrics']))
            if not hsm:
                hsm = exome_qc_files['HsMetrics']
            else:
                hsm = hsm[0]
            if os.path.isfile(hsm):
                header_add_fields = hs_metrics(hsm)
                metrics_header = met_wgs_header + header_add_fields
                # metrics_header = list(results.keys())
                # metrics_header.sort()
            else:
                metrics_header = met_wgs_header

            os.chdir(working_dir)
            write_results(results, cwd_metrics_outfile, metrics_header)

            # process individual libraries
            if args.l:

                sample_library = []
                os.chdir(info[3] + '/results')

                mdm = glob.glob('*{}*'.format(exome_qc_files['mark_dups_metrics']))
                if not mdm:
                    sys.exit('No mark_dups_metrics found')

                with open(mdm[0], 'r') as getmettxt:
                    infile_reader = csv.reader(getmettxt, delimiter='\t')
                    for line in infile_reader:
                        if line and 'lib' in line[0]:
                            sample_library.append(line[0])

                print()
                print(info[1], info[4], 'Libraries:')
                for sample in sample_library:
                    print(sample)
                    os.chdir(info[3] + '/results')
                    results['Library'] = sample

                    # update mark_dup_metrics for individual library
                    with open(mdm[0], 'r') as getmettxt:
                        infile_reader = csv.reader(getmettxt, delimiter='\t')
                        for line in infile_reader:
                            if line and sample in line[0]:
                                results['PERCENT_DUPLICATION'] = line[8]

                    # update alignment summary metrics for individual library
                    asm = glob.glob('*{}*'.format(exome_qc_files['AlignmentSummaryMetrics']))
                    if not asm:
                        sys.exit('No AlignmentSummaryMetrics file found.')
                    with open(asm[0], 'r') as aligntxt:
                        infile_reader = csv.reader(aligntxt, delimiter='\t')
                        pf_aligned_bases = float()
                        for line in infile_reader:
                            if sample in line:
                                if 'FIRST_OF_PAIR' in line:
                                    results['FOP: PF_MISMATCH_RATE'] = line[12]
                                if 'SECOND_OF_PAIR' in line:
                                    results['SOP: PF_MISMATCH_RATE'] = line[12]
                                if 'PAIR' in line and not '_' in line:
                                    pf_aligned_bases = int(line[7])
                                    results['TOTAL_READS'] = line[1]
                                    results['PF_READS'] = line[2]
                                    results['PF_READS_ALIGNED'] = line[5]
                                    results['PF_ALIGNED_BASES'] = line[7]
                                    results['PF_HQ_ALIGNED_Q20_BASE'] = line[10]
                                    results['PCT_ADAPTER'] = line[23]

                    # update HsMetrics.txt for individual library
                    hsm = glob.glob('*{}*'.format(exome_qc_files['HsMetrics']))
                    if not hsm:
                        sys.exit('No HsMetrics file found.')
                    with open(hsm[0], 'r') as hasmettxt:
                        infile_reader = csv.reader(hasmettxt, delimiter='\t')
                        for line in infile_reader:
                            if 'BAIT_SET' in line:
                                hs_metrics_header = line
                            if sample in line:
                                hs_metrics_data_two = line
                                hs_metrics_dict = dict(zip(hs_metrics_header, hs_metrics_data_two))
                        for metric in hs_metrics_dict:
                            results[metric] = hs_metrics_dict[metric]

                    # update InsertSizeMetrics.txt for individual library
                    ism = glob.glob('*{}*'.format(exome_qc_files['InsertSizeMetrics']))
                    if not ism:
                        sys.exit('No InsertSizeMetrics file found.')
                    with open(ism[0], 'r') as inserttxt:
                        infile_reader = csv.reader(inserttxt, delimiter='\t')
                        for line in infile_reader:
                            if 'MEAN_INSERT_SIZE' in line:
                                mean_header_position = line.index('MEAN_INSERT_SIZE')
                            if sample in line:
                                results['MEAN_INSERT_SIZE'] = line[mean_header_position]
                                results['STANDARD_DEVIATION'] = line[mean_header_position+1]

                    # write library results to outfile
                    os.chdir(working_dir)
                    write_library_results(results, cwd_metrics_library_outfile, metrics_header)

        else:
            print('{} {} Build Failed'.format(info[0], info[1]))

if args.e:
    print('----------')
    print('cwl.metrics.parser.py generation complete.\nRunning exome_report on cwl exome metrics files.')
    print('----------')
    subprocess.run(["/gscuser/ltrani/Desktop/python/bin/python3", "/gscmnt/gc2783/qc/bin/aw/exome_report.py"])

if args.wgs:
    print('----------')
    print('cwl.metrics.parser.py generation complete.\nRunning wgs_report on cwl wgs metrics files.')
    print('----------')
    subprocess.run(["/gscuser/ltrani/Desktop/python/bin/python3", "/gscmnt/gc2783/qc/bin//aw/wgs_report.py"])
