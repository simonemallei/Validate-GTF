# Validate GTF
# Student Code Number : 844659
# Name : Simone
# Surname : Mallei

import re
from datetime import datetime

# Given a string containing the ninth field of a gtf record, it returns
# the transcript_id and the gene_id
def get_transcript_and_gene(attrs):
    transcript_id = re.findall('transcript_id\s+"([^"]*)";', attrs)[0]
    gene_id = re.findall('gene_id\s+"([^"]*)";', attrs)[0]
    return (transcript_id, gene_id)

# Given a string containing a gtf record, its index in the file and
# the dictionary containing all the errors that has been found, it verifies 
# intra-row constraints and modify the dictionary depending on the 
# violations of gtf_row
def validate_row(gtf_row, index_row, errors):
    num_fields = 9
        
    # It removes the comment from the row
    # Checks number of fields
    gtf_row = re.findall('([^#]*)#?', gtf_row)[0] 
    gtf_fields = gtf_row.split('\t')
        
    if len(gtf_fields) != num_fields:
        errors['num_fields'].append(index_row)
        return

    # Verifies <feature> constraints
    feature_list = ['CDS', 'start_codon', 'stop_codon', '5UTR', '3UTR', 'inter', 'inter_CNS', 'intron_CNS', 'exon']
    gtf_feature = gtf_fields[2]
        
    # Verifies if the feature is to ignore or not
    if gtf_feature not in feature_list:
        errors['feature_name'].append(index_row)
        return
        
    # Verifies <start> and <end> constraints
    try:
        start_coord = int(gtf_fields[3])
        if start_coord < 1:
            error['start'].append(index_row)
    except:
        errors['start'].append(index_row)
    if index_row not in (errors['start']):
        try:
            end_coord = int(gtf_fields[4])
            if end_coord < 1:
                errors['end'].append(index_row)
            elif start_coord > end_coord:
                errors['start_end'].append(index_row)
            elif gtf_feature == 'start_codon' and end_coord - start_coord > 2:
                errors['start_codon_len'].append(index_row)
            elif gtf_feature == 'stop_codon' and end_coord - start_coord > 2:
                errors['stop_codon_len'].append(index_row)
        except:
            errors['end'].append(index_row)

    # Verifies <score> constraints
    gtf_score = gtf_fields[5]
    if gtf_score != '.':
        try :
            float(gtf_score)
        except :
            errors['score'].append(index_row)

    # Verifies <strand> constraints
    gtf_strand = gtf_fields[6]
    if gtf_strand not in ('+', '-'):
        errors['strand'].append(index_row)
            
    # Verifies <frame> constraints
    if gtf_feature in ('CDS', 'start_codon', 'stop_codon'):
        try :
            gtf_frame = int(gtf_fields[7])
            if gtf_frame < 0 or gtf_frame > 2:
                errors['frame'].append(index_row)
        except :
            errors['frame'].append(index_row)
    elif gtf_fields[7] != '.':
        errors['frame'].append(index_row)

    # Verifies <attributes> constraints
    attrs = gtf_fields[8].rstrip()
    try:
        transcript_id = re.findall('transcript_id\s"([^"]*)";', attrs)[0]
        gene_id = re.findall('gene_id\s"([^"]*)";', attrs)[0]
        attributes_start = 'gene_id "'+gene_id+'"; '
        attributes_start += 'transcript_id "'+transcript_id+'";'
        if not(gtf_fields[8].startswith(attributes_start)):
            errors['attributes'].append(index_row)
            return
        pat_num = '(([-,+]?\d+(.\d*)([e,E][-,+]?\d+)?)|[N,n][A,a][N,n])'
        pat_text = '"[^"]*"'
        pat_attr = '[^"\s]+\s('+pat_num+'|'+pat_text+');'
        pattern = '(('+pat_attr+'\s)*('+pat_attr+'))'
        attributes = re.findall(pattern, attrs)[0][0]

        # Verifies that if <feature> is == 'inter' or 'inter_CNS', then
        # transcript_id has to be empty
        if gtf_feature in ('inter', 'inter_CNS'):
            if transcript_id != '':
                errors['inter_transcript'].append(index_row+1)
        elif attributes != attrs or gene_id == '' or transcript_id == '':
            errors['attributes'].append(index_row)
    except:
        errors['attributes'].append(index_row)

# Given an errors' dictionary and the name of the file that should have
# these errors, reports in a file which violations occurred in the file.
def write_report(errors, file_name):
    error_string = {
        'num_fields': 'FieldsError: Wrong number of fields in row: ',
        'start': 'StartError: Start value is less than 1 in row: ',
        'end': 'EndError: End value is less than 1 in row: ',
        'start_end': 'StartEndError: Start value greater than end value in row: ',
        'start_codon_len': 'StartCodonError: Wrong length of start codon in row: ',
        'stop_codon_len': 'StopCodonError: Wrong length of stop codon in row: ',
        'feature_name': 'FeatureError: Wrong feature name in row: ',
        'inter_transcript': 'TranscriptError: Transcript not empty in "inter" or "inter_CNS" feature in row: ',
        'score': 'ScoreError: Score not integer, float or "." in row: ',
        'strand': 'StrandError: Strand is not "+" or "-" in row: ',
        'frame': 'FrameError: Frame has value not available in row: ',
        'attributes': 'AttributesError: Wrong attributes in row: ',
        'strand_interrow': 'StrandError: Both strands used in: ',
        'missing_CDS': 'CDSError: Missing CDS feature in: ',
        'missing_start_codon': 'StartCodonError: Missing start_codon feature in: ',
        'missing_stop_codon': 'StopCodonError: Missing stop_codon feature in: ',
        'CDS_interrow': 'CDSError: CDS intervals set wrong in : ',
        'start_codon_interrow': 'StartCodonError: start_codon intervals set wrong in: ',
        'stop_codon_interrow': 'StopCodonError: stop_codon intervals set wrong in: ',
        'CDS_start_codon': 'CDSStartCodonError: CDS does not start with start_codon in: ',
        'CDS_stop_codon': 'CDSStopCodonError: stop_codon is not after CDS in: '}
    dateString = str(datetime.today())
    file_no_ext = re.findall('(\w*).gtf', file_name)[0]
    output_file = open('reports/report-'+file_no_ext+'.txt', 'w')
    output_file.write('REPORT OF FILE: '+file_name+'\nTime of report: '+dateString+'\n\n')
    # Prints violations in stdout and it writes them in the output file
    for error_type in errors:
        for element in errors[error_type]:
            print(error_string[error_type]+str(element))
            output_file.write(error_string[error_type]+str(element)+'\n')
    output_file.close()

# Given a path of the directory containing the file to validate and the file's name
# it verifies constraints in the entire .gtf file (intra-row and inter-row constraints)
def validate_file(gtf_dir, gtf_file_name):
    
    with open(gtf_dir+gtf_file_name, 'r') as gtf_input_file:
        gtf_file_rows = gtf_input_file.readlines()

    errors = {'num_fields': [],
              'start': [],
              'end': [],
              'start_end': [],
              'start_codon_len': [],
              'stop_codon_len': [],
              'feature_name': [],
              'inter_transcript': [],
              'score': [],
              'strand': [],
              'frame': [],
              'attributes': [],
              'strand_interrow': set(),
              'missing_CDS': set(),
              'missing_start_codon': set(),
              'missing_stop_codon': set(),
              'CDS_interrow': set(),
              'start_codon_interrow': set(),
              'stop_codon_interrow': set(),
              'CDS_start_codon': set(),
              'CDS_stop_codon': set()}
    

    indexes = range(len(gtf_file_rows))
    # Validate each row (verifying only intra-row violations)
    for (gtf_row, index_row) in zip(gtf_file_rows, indexes):
        validate_row(gtf_row, index_row+1, errors)

    strand_dict = {}
    features_dict = {}

    # Verifies inter-row violations
    valid_rows = []
    for index_row in indexes:
        is_valid_row = True
        for error_type in errors:
            if index_row+1 in errors[error_type]:
                is_valid_row = False
        if is_valid_row:
            valid_rows.append(index_row)
            
    for index_row in valid_rows:
        gtf_row = re.findall('([^#]*)#?', gtf_file_rows[index_row])[0] 
        gtf_fields = gtf_row.split('\t')
        gtf_feature = gtf_fields[2]
        gtf_strand = gtf_fields[6]
        transcript_id, gene_id = get_transcript_and_gene(gtf_fields[8])
        
        # Maps 'CDS', 'start_codon', 'stop_codon' features to lists depending
        # on their transcript_id
        
        if gtf_feature in ('CDS', 'start_codon', 'stop_codon'):
            gtf_start = int(gtf_fields[3])
            gtf_end = int(gtf_fields[4])
            gtf_frame = int(gtf_fields[7])
            features_dict[gene_id] = features_dict.get(gene_id, {})
            features_dict[gene_id][transcript_id] = features_dict[gene_id].get(transcript_id, {})
            features_dict[gene_id][transcript_id][gtf_feature] = features_dict[gene_id][transcript_id].get(gtf_feature, [])
            features_dict[gene_id][transcript_id][gtf_feature].append((gtf_start, gtf_end, gtf_frame))

        # Verifies <strand> constraints for each gene_id
        if strand_dict.get(gene_id, '') == '':
            strand_dict[gene_id] = gtf_strand
        elif strand_dict[gene_id] != gtf_strand:
            errors['strand_interrow'].add(gene_id)
            
    required_features = ('CDS', 'start_codon', 'stop_codon')
    for gene_id in features_dict.keys():
        if gene_id not in errors['strand_interrow']:
            for transcript_id in features_dict[gene_id].keys():
                
                # Violations to verify in each transcript
                # Sort features depending on <strand>

                alreadyViolated = False

                for gtf_feature in required_features:
                    curr_feature_list = features_dict[gene_id][transcript_id].get(gtf_feature, [])
                    sorted_list = sorted(curr_feature_list)
                    if strand_dict[gene_id] == '-':
                        sorted_list.reverse()
                    features_dict[gene_id][transcript_id][gtf_feature] = sorted_list
                    
                    if sorted_list == [] and gtf_feature in required_features:
                        errors['missing_'+gtf_feature].add((gene_id, transcript_id))
                        alreadyViolated = True
                    
                if alreadyViolated:
                    continue;
                
                # Verifies that every 'CDS' feature does not intersect with
                # other 'CDS' feature, that every frame annotated is correct and
                # the total length is correct.
                
                for gtf_feature in ('CDS', 'start_codon', 'stop_codon'):
                    
                    codon_nucleotides = features_dict[gene_id][transcript_id][gtf_feature]
                    codon_len = 0
                    last_nucleotide = 0
                    frame_calc = 0
                    for (curr_start, curr_end, curr_frame) in codon_nucleotides:
                        if curr_frame != frame_calc:
                            errors[gtf_feature + '_interrow'].add((gene_id, transcript_id))
                        elif (gene_id, transcript_id) not in errors[gtf_feature + '_interrow']:
                            if codon_len == 0:
                                if strand_dict[gene_id] == '+':
                                    codon_len += curr_end - curr_start + 1
                                    last_nucleotide = curr_end
                                else:
                                    codon_len += curr_end - curr_start + 1
                                    last_nucleotide = curr_start
                            elif (gene_id, transcript_id) not in errors[gtf_feature + '_interrow']:
                                if strand_dict[gene_id] == '+':
                                    if curr_start <= last_nucleotide:
                                        errors[gtf_feature + '_interrow'].add((gene_id, transcript_id))
                                    else:
                                        codon_len += curr_end - curr_start + 1
                                        last_nucleotide = curr_end
                                else:
                                    if curr_end >= last_nucleotide:
                                        errors[gtf_feature + '_interrow'].add((gene_id, transcript_id))
                                    else:
                                        codon_len += curr_end - curr_start + 1
                                        last_nucleotide = curr_start
                            frame_calc = (3 - ((curr_end - curr_start + 1 - frame_calc) % 3)) % 3
                    if (gene_id, transcript_id) not in errors[gtf_feature + '_interrow']:
                        if gtf_feature != 'CDS' and codon_len != 3:
                            errors[gtf_feature + '_interrow'].add((gene_id, transcript_id))
                        elif gtf_feature == 'CDS' and codon_len % 3 != 0:
                            errors[gtf_feature + '_interrow'].add((gene_id, transcript_id))


                # Verifies that CDS starts with start_codon
                if (gene_id, transcript_id) not in errors['CDS_interrow']:
                    codon_intervals = len(features_dict[gene_id][transcript_id]['start_codon'])
                    CDS_list = features_dict[gene_id][transcript_id]['CDS']
                    if (gene_id, transcript_id) not in errors['start_codon_interrow']:
                        start_codon_list = features_dict[gene_id][transcript_id]['start_codon']
                        if len(start_codon_list) > len(CDS_list):
                            errors['CDS_start_codon'].add((gene_id, transcript_id))
                        for index in range(codon_intervals):
                            if index == codon_intervals -1:
                                if strand_dict[gene_id] == '+':
                                    if start_codon_list[index][0] != CDS_list[index][0] or start_codon_list[index][1] > CDS_list[index][1]:
                                        errors['CDS_start_codon'].add((gene_id, transcript_id))
                                elif strand_dict[gene_id] == '-':
                                    if start_codon_list[index][1] != CDS_list[index][1] or start_codon_list[index][0] < CDS_list[index][0]:
                                        errors['CDS_start_codon'].add((gene_id, transcript_id))
                            elif start_codon_list[index] != CDS_list[index]:
                                errors['CDS_start_codon'].add((gene_id, transcript_id))

                    # Verifies that stop_codon's positions are after CDS    
                    if (gene_id, transcript_id) not in errors['stop_codon_interrow']:
                        stop_codon_list = features_dict[gene_id][transcript_id]['stop_codon']
                        last_CDS = CDS_list[-1]
                        if strand_dict[gene_id] == '+':
                            if stop_codon_list[0][0] <= CDS_list[-1][1]:
                                errors['CDS_stop_codon'].add((gene_id, transcript_id))
                        else:
                            if stop_codon_list[0][1] >= CDS_list[-1][0]:
                                errors['CDS_stop_codon'].add((gene_id, transcript_id))

    write_report(errors, gtf_file_name)
