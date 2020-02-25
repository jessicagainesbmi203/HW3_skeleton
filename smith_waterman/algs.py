import numpy as np
import pandas as pd
from .io import read_sub_matrix, read_pairs, read_sequence
import re
import matplotlib.pyplot as plt

def smithwaterman(seq1,seq2,sub_matrix_file,gap_open,gap_extend):
    """
    Perform sequence alignment according to the Smith-Waterman Algorithm
    Inputs: seq1 and seq2: two sequence objects containing sequences to compare
        sub_matrix_file: filepath to desired substitution matrix
        gap_open: penalty subtracted when a new gap is opened
        gap_extend: penalty subtracted when a gap is extended
    Outputs: a dictionary defining three outputs for easy identification
        aligned_seq1 : the segment of seq1, with inserted gaps, that had optimal local alignment with seq2
        aligned_seq2 : the segment of seq2, with inserted gaps, that had optimal local alignment with seq1
        alignment_score : the sum of all elements of the scoring matrix along the alignment path
    """
    # check if either sequence is empty
    if seq1.sequence == "" or seq2.sequence == "": 
        return {'aligned_seq1':'', 'aligned_seq2':'', 'alignment_score':0}
    scoring_matrix_filepath = 'scoring_matrix/' + seq1.name + '_' + seq2.name + '_' + sub_matrix_file + '_' + str(gap_open) + '_' + str(gap_extend) + '.csv'
    # save and read in scoring matrix to save time
    # if time before assignment is due, work on efficiency in this section of the code
    try:
        scoring_matrix = pd.read_csv(scoring_matrix_filepath,header=0,index_col=0)
        scoring_matrix.columns = range(len(scoring_matrix.shape[1]))
    except:
        print('create scoring matrix')
        scoring_matrix = create_scoring_matrix(seq1,seq2,sub_matrix_file,gap_open,gap_extend)
        scoring_matrix.to_csv(scoring_matrix_filepath)
    alignment = find_alignment(seq1,seq2,scoring_matrix)
    return alignment

def create_scoring_matrix(seq1,seq2,sub_matrix_file,gap_open,gap_extend):
    """
    Create a scoring matrix to determine optimal alignment of the two sequences, based on probability of substitution and gap penalties
    Inputs: seq1 and seq2: two sequence objects containing sequences to compare
        sub_matrix_file: filepath to desired substitution matrix
        gap_open: penalty subtracted when a new gap is opened
        gap_extend: penalty subtracted when a gap is extended
    Output: scoring matrix to trace back to optimal sequence alignment
    """
    sub_matrix = read_sub_matrix(sub_matrix_file)
    # initialize scoring matrix, set first row and column to 0
    scoring_matrix = pd.DataFrame(index=range(len(seq1.sequence)+1),columns=range(len(seq2.sequence)+1))
    scoring_matrix.iloc[:,0] = 0
    scoring_matrix.iloc[0,:] = 0
    # fill in scoring matrix with max score of gap in seq1, gap in seq2, no gap, or 0
    k = 0
    l = 0
    for i in range(1,len(scoring_matrix.index),1):
        res1 = seq1.sequence[i-1]
        for j in range(1,len(scoring_matrix.columns),1):
            res2 = seq2.sequence[j-1]
            # find score if no gaps -- the score 
            no_gaps = int(scoring_matrix.loc[i-1,j-1]) + int(sub_matrix.loc[res1,res2])
            # comment out this piece of code -- it takes into account all possible gap lengths but it takes too long
                #gap_seq1 = 0
                #for k in range(1,i,1):
                #    score = scoring_matrix.loc[i-k,j] - gap_open - (k-1)*gap_extend
                #    if (score > gap_seq1):
                #        gap_seq1 = score
                #gap_seq2 = 0
                #for l in range(1,j,1):
                #    score = scoring_matrix.loc[i,j-l] - gap_open - (l-1)*gap_extend
                #    if (score > gap_seq2):
                #        gap_seq2 = score
            # instead predict how long the gap will be -- keep track of what the traceback will choose
            gap_seq1 = scoring_matrix.loc[i-1,j] - gap_open - k*gap_extend
            gap_seq2 = scoring_matrix.loc[i,j-1] - gap_open - l*gap_extend
            scoring_matrix.loc[i,j] = max(no_gaps,gap_seq1,gap_seq2,0)
            if max(no_gaps,gap_seq1,gap_seq2,0) == gap_seq1:
                k += 1
                l = 0
            elif max(no_gaps,gap_seq1,gap_seq2,0) == gap_seq2:
                l += 1
                k = 0
            else:
                k = 0
                l = 0
    return scoring_matrix

def find_alignment(seq1,seq2,scoring_matrix):
    """
    Using the previously generated scoring matrix, retrace through the scoring matrix to find the best local alignment.
    Start at the highest score and trace up and left until a score of zero is reached.
    Inputs: seq1 and seq2: two sequence objects containing sequences to compare
        scoring_matrix: the output of create_scoring_matrix() which stores scores related to the probability of the sequences being related
    Outputs: a dictionary defining three outputs for easy identification
        aligned_seq1 : the segment of seq1, with inserted gaps, that had optimal local alignment with seq2
        aligned_seq2 : the segment of seq2, with inserted gaps, that had optimal local alignment with seq1
        alignment_score : the sum of all elements of the scoring matrix along the alignment path
    """
    # Find the highest score in the scoring matrix, use this location as a starting point
    row = scoring_matrix.max(axis=1).idxmax()
    col = scoring_matrix.max(axis=0).idxmax()
    # Follow the path of highest scores toward the beginning of the two sequences
    # Track sequence alignment and max score.
    aligned_seq1 = list()
    aligned_seq2 = list()
    max_score = scoring_matrix.loc[row,col]
    aligned_seq1.append(seq1.sequence[row-1])
    aligned_seq2.append(seq2.sequence[col-1])
    # continue until zero score is found
    while (scoring_matrix.loc[row,col] > 0):
        no_gap = scoring_matrix.loc[row-1,col-1]
        gap_seq1 = scoring_matrix.loc[row,col-1]
        gap_seq2 = scoring_matrix.loc[row-1,col]
        # If adding no gap is optimal
        if no_gap >= gap_seq1 and no_gap >= gap_seq2 :
            row = row-1
            col = col-1
            if scoring_matrix.loc[row,col] > 0:
                aligned_seq1.append(seq1.sequence[row-1])
                aligned_seq2.append(seq2.sequence[col-1])
        # Add gap in sequence 1 is optimal
        if (gap_seq1 > no_gap and gap_seq1 >= gap_seq2):
            col = col-1
            if scoring_matrix.loc[row,col] > 0:
                aligned_seq1.append('-')
                aligned_seq2.append(seq2.sequence[col-1])
        # Add gap in sequence 2 is optimal
        if gap_seq2 > no_gap and gap_seq2 > gap_seq1:
            row = row-1
            if scoring_matrix.loc[row,col] > 0:
                aligned_seq1.append(seq1.sequence[row-1])
                aligned_seq2.append('-')
    # reverse sequence order to match input -- it was built backwards
    aligned_seq1.reverse()
    aligned_seq2.reverse()
    return {'aligned_seq1':aligned_seq1, 'aligned_seq2':aligned_seq2, 'alignment_score':max_score}

def score(aligned_seq1,aligned_seq2,sub_matrix_file,gap_open,gap_extend):
    """
    Take into account substitution matrix scores and gap penalties to find the alignment score of an alignment
    Inputs: aligned_seq1, aligned_seq2: sequences with gaps included to create optimal alignment between the two sequences
        sub_matrix_file: filepath to desired substitution matrix
        gap_open: penalty subtracted when a new gap is opened
        gap_extend: penalty subtracted when a gap is extended
    Outputs: alignment score of the two sequences
    """
    sub_matrix = read_sub_matrix(sub_matrix_file)
    cumulative_score = 0
    length = min(len(aligned_seq1),len(aligned_seq2))
    gap_open_flag_seq1 = 0
    gap_open_flag_seq2 = 0
    sub_score = 0
    # iterate through each pair of sequences
    for i in range(length):
        res1 = aligned_seq1[i]
        res2 = aligned_seq2[i]
        # find gap penalty for a gap in sequence 1
        if res1 == '-':
            # if this is the first gap character, use the gap opening penalty
            if not gap_open_flag_seq1:
                sub_score -= gap_open
                gap_open_flag_seq1 = 1
            # if this is not the first gap character, use the gap extension penalty
            else:
                sub_score -= gap_extend
            gap_open_flag_seq2 = 0
        if res2 == '-':
            if not gap_open_flag_seq2:
                sub_score -= gap_open
                gap_open_flag_seq2 = 1
            else:
                sub_score -= gap_extend
            gap_open_flag_seq1 = 0
        # if there are no gaps recorded, reset gaps and find substitution matrix score
        if res1 != '-' and res2 != '-':
            gap_open_flag_seq1 = 0
            gap_open_flag_seq2 = 0
            sub_score += sub_matrix.loc[res1,res2]
            # smith-waterman does not allow negative scores
            cumulative_score += max(sub_score,0)
            sub_score = 0
    return cumulative_score

def get_false_pos_rate(sub_matrix_file,gap_open,gap_extend,true_filepath, false_filepath, true_pos_rate, normalize=False):
    """
    Use the lists of positive and negative pairs to evaluate the algorithm for each substitution matrix and set of gap penalties
    Inputs: sub_matrix_file: filepath to desired substitution matrix
        gap_open: penalty subtracted when a new gap is opened
        gap_extend: penalty subtracted when a gap is extended
        true_filepath: filepath to the document containing a list of pairs of related sequences
        false_filepath: filepath to the document containing a list of pairs of unrelated sequences
    Outputs: false_pos_rate: the proportion of pairs of unrelated proteins that are flagged as related at the alignment score 
        that creates a 70% detection of true positives
    """
    # to save time, read in the dataframe of positive pairs and their scores if it exists
    true_df_filepath = 'Pospairs_scores' + '_' + sub_matrix_file + '_' + str(gap_open) + '_' + str(gap_extend) + '.csv'
    try:
        true_df = pd.read_csv(true_df_filepath)
    except:
        true_df = pairs_sorted_dataframe(true_filepath, sub_matrix_file, gap_open, gap_extend, true_df_filepath)
    # to save time, read in the dataframe of positive pairs and their scores if it exists
    false_df_filepath = 'Negpairs_scores' + '_' + sub_matrix_file + '_' + str(gap_open) + '_' + str(gap_extend) + '.csv'
    try:
        false_df = pd.read_csv(false_df_filepath)
    except:
        false_df = pairs_sorted_dataframe(false_filepath, sub_matrix_file, gap_open, gap_extend, false_df_filepath)
    # find the threshold at 30% * data length from the lowest score in the array sorted ascendiing
    if normalize:
        score_col = 'normalized_score'
    else:
        score_col = 'score'
    true_df = true_df.sort_values(by=score_col)
    false_df = false_df.sort_values(by=score_col)
    threshold = true_df.loc[np.floor(true_df.shape[0] * (1-true_pos_rate)),score_col]
    print(threshold)
    # apply threshold to negative pairs and find the number above threshold
    false_pos = false_df[false_df[score_col] > threshold]
    print(false_pos)
    # find the proportion of negative pairs above threshold
    false_pos_rate = false_pos.shape[0] / false_df.shape[0]
    return false_pos_rate
    
def get_true_pos_rate(sub_matrix_file,gap_open,gap_extend,true_filepath, false_filepath, false_pos_rate, normalize=False):
    """
    Use the lists of positive and negative pairs to evaluate the algorithm for each substitution matrix and set of gap penalties
    Inputs: sub_matrix_file: filepath to desired substitution matrix
        gap_open: penalty subtracted when a new gap is opened
        gap_extend: penalty subtracted when a gap is extended
        true_filepath: filepath to the document containing a list of pairs of related sequences
        false_filepath: filepath to the document containing a list of pairs of unrelated sequences
    Outputs: true_pos_rate: the proportion of pairs of related proteins that are flagged as related at the alignment score 
        that creates a given detection of false positives
    PS Sorry for the copy-pasted code, with more time I would somehow merge this with get_false_pos_rate
    """
    # to save time, read in the dataframe of positive pairs and their scores if it exists
    true_df_filepath = 'Pospairs_scores' + '_' + sub_matrix_file + '_' + str(gap_open) + '_' + str(gap_extend) + '.csv'
    try:
        true_df = pd.read_csv(true_df_filepath)
    except:
        true_df = pairs_sorted_dataframe(true_filepath, sub_matrix_file, gap_open, gap_extend, true_df_filepath)
    # to save time, read in the dataframe of positive pairs and their scores if it exists
    false_df_filepath = 'Negpairs_scores' + '_' + sub_matrix_file + '_' + str(gap_open) + '_' + str(gap_extend) + '.csv'
    try:
        false_df = pd.read_csv(false_df_filepath)
    except:
        false_df = pairs_sorted_dataframe(false_filepath, sub_matrix_file, gap_open, gap_extend, false_df_filepath)
    # find the threshold at 30% * data length from the lowest score in the array sorted ascendiing
    if normalize:
        score_col = 'normalized_score'
    else:
        score_col = 'score'
    true_df = true_df.sort_values(by=score_col)
    print(true_df)
    false_df = false_df.sort_values(by=score_col)
    print(false_df)
    threshold = false_df.loc[np.floor(false_df.shape[0] * (1-false_pos_rate)),score_col]
    # apply threshold to negative pairs and find the number above threshold
    true_pos = true_df[true_df[score_col] > threshold]
    # find the proportion of negative pairs above threshold
    true_pos_rate = true_pos.shape[0] / true_df.shape[0]
    return true_pos_rate
    
def pairs_sorted_dataframe(pairs_filepath, sub_matrix_file, gap_open, gap_extend, output_filepath):
    """
    For a list of pairs of sequences, find the alignment score for each pair and store in a dataframe, with pairs sorted by alignment score
    Inputs: pairs_filepath: filepath to the document containing the pairs
        sub_matrix_file: filepath to desired substitution matrix
        gap_open: penalty subtracted when a new gap is opened
        gap_extend: penalty subtracted when a gap is extended
        output_file_path: filepath to save sorted dataframe
    Outputs: dataframe with pairs of sequences sorted by alignment score
    """
    pairs = read_pairs(pairs_filepath)
    df = pd.DataFrame(columns=('seq1','seq2','score'),index=range(len(pairs)))
    i = 0
    for pair in pairs:
        print(pair)
        alignment = smithwaterman(pair[0],pair[1],sub_matrix_file,gap_open,gap_extend)
        print(alignment)
        alignment_score = score(alignment.get('aligned_seq1'), alignment.get('aligned_seq2'), sub_matrix_file,gap_open,gap_extend)
        df.loc[i,'seq1'] = pair[0]
        df.loc[i,'seq2'] = pair[1]
        df.loc[i,'score'] = alignment_score
        length = min(len(pair[0].sequence),len(pair[1].sequence))
        df.loc[i,'normalized_score'] = alignment_score / length
        i += 1
        print(alignment_score)
    print(df)
    df = df.sort_values(by='score')
    df.to_csv(output_filepath)
    print(df)
    return df
    
def roc(sub_matrix_file,gap_open,gap_extend,true_filepath, false_filepath, false_pos_rate_list, normalize):
    """
    Create a plot that compares true positive rate to false positive rate
    Inputs: sub_matrix_file, a dataframe with substitution scores for each pair of residues
        gap_open: penalty subtracted when a new gap is opened
        gap_extend: penalty subtracted when a gap is extended
        true_filepath: filepath to the document containing a list of pairs of related sequences
        false_filepath: filepath to the document containing a list of pairs of unrelated sequences
        true_pos_rate_list: list a values of true positive rates to find corresponding false positive rates
        normalize: boolean, whether or not the score is divided by the length of the shorter sequence
    """
    y = list()
    for false_pos_rate in false_pos_rate_list:
        y.append(get_true_pos_rate('BLOSUM50',10,2,'Pospairs.txt','Negpairs.txt',false_pos_rate, normalize))
    plt.figure(figsize=(5,5))
    plt.plot(false_pos_rate_list,y)
    plt.xlabel('fraction false positives')
    plt.ylabel('fraction true positives')
    plt.title('ROC curve')
    plt.xlim([0,1])
    plt.ylim([0,1])
    plt.show()
    return y

'''
Part II
'''

def optimize(sub_matrix_file, true_filepath, false_filepath, fp_rates):
    population = mutate(sub_matrix_file,100)
    df = pd.DataFrame(columns=('score','matrix'))
    i = 0
    for matrix in population:
        score = objective_function(fp_rates, sub_matrix_file)
        store.loc[i,'score'] = score
        store.loc[i,'matrix'] = matrix
        i += 1
    df = df.sort_values(by='score', ascending=False)
    
def objective_function(fp_rates,sub_matrix_file):
    sum = 0
    for rate in fp_rates:
        tp_rate = get_true_pos_rate(sub_matrix_file,10,2,'Pospairs.txt','Negpairs.txt',rate, normalize=False)
        sum += tp_rate
    return sum
        
def mutate(sub_matrix_file,n):
    sub_matrix = read_sub_matrix(sub_matrix_file)
    population = list()
    for i in n:
        new_matrix = sub_matrix.copy()
        for j in range(new_matrix.shape[0]):
            for k in range(j):
                new_element = new_matrix.iloc[j,k] + (0.5*np.random.random() - 0.25)
                new_matrix.iloc[j,k] = new_element
                new_matrix.iloc[k,j] = new_element
        population.append(new_matrix)
    return population
        
    