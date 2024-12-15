import streamlit as st
from PIL import Image
import numpy as np
import pandas as pd
from funcations import *

if "page" not in st.session_state:
    st.session_state.page = "Home"

def parse_fasta(fasta_file):
    if fasta_file is None:
        return None
    sequence = []
    for line in fasta_file:
        line = line.decode("utf-8").strip()
        if not line.startswith(">"):
            sequence.append(line)
    return "".join(sequence)
    
st.set_page_config(page_title="Bioinformatics", layout="wide", initial_sidebar_state="collapsed")

with st.sidebar:
    if st.button('Home'):
        st.session_state.page = "Home"

    if st.button("Complement"):
        st.session_state.page = "Complement"

    if st.button("Reverse"):
        st.session_state.page = "Reverse"

    if st.button('Reverse Complement'):
        st.session_state.page = 'Reverse_Complement'

    if st.button('Translation Table'):
        st.session_state.page = 'Translation_Table'

    if st.button('Match'):
        st.session_state.page = 'Match'

    if st.button('Bad char'):
        st.session_state.page = 'Bad'

    if st.button('Index Sorted'):
        st.session_state.page = 'Index'

    if st.button('Query'):
        st.session_state.page = 'Query'

    if st.button('Suffix array construction'):
        st.session_state.page = 'Suffix_array_construction'

    if st.button('Overlap'):
        st.session_state.page = 'Overlap'

    if st.button('Native overlap'):
        st.session_state.page = 'native_overlap'

if st.session_state.page == "Home":
    st.title("Introduction to Bioinformatics")
    st.write('A DNA sequence is made up of four nucleotides: A (Adenine), T (Thymine), C (Cytosine), and G (Guanine).')
    st.write('**Select what you want from side bar**')

if st.session_state.page == "Complement":
    st.title('**Complement of a DNA Sequence**')  
    st.write('The complement of a DNA strand is formed by replacing each nucleotide with its complementary base:')
    s = """
    * A pairs with T.
    * T pairs with A.
    * C pairs with G.
    * G pairs with C.
    """
    st.write(s)
    seq = st.text_area("Enter Sequence for DNA:", placeholder="Type here...")  
    if st.button("Submit"):
        if seq:
            st.write(f'Complement Sequence is: ', Complement(seq))
    fasta_file = st.file_uploader("Upload a FASTA or FNA file", type=["fasta", "fna"])
    if fasta_file:
        sequence = parse_fasta(fasta_file)
        if sequence:
            st.write(f'Sequence is: {sequence}')
            # st.write("### Complement Sequence:")
            st.write(f'Complement Sequence is: {Complement(sequence)}')

if st.session_state.page == "Reverse":
    st.title('**Reverse of a DNA Sequence**')
    st.write('The reverse of a DNA sequence is simply the original sequence written in reverse order,\
              without changing any of the nucleotide bases.')
    seq = st.text_area("Enter Sequence for DNA:", placeholder="Type here...")
    if st.button("Submit"):
        if seq:
            st.write(f'Reverse Sequence is: ', Reverse(seq))
    fasta_file = st.file_uploader("Upload a FASTA or FNA file", type=["fasta", "fna"])
    if fasta_file:
        sequence = parse_fasta(fasta_file)
        if sequence:
            st.write(f'Sequence is: {sequence}')
            # st.write("### Reverse Sequence:")
            st.write(f'Reverse Sequence is: {Reverse(sequence)}')

if st.session_state.page == 'Reverse_Complement':
    st.title('**Reverse Complement of a DNA Sequence**')
    st.write('First finds the complement, then reverses the sequence.')
    seq = st.text_area("Enter Sequence for DNA:", placeholder="Type here...")
    if st.button("Submit"):
        if seq:
            st.write(f'Reverse Complement Sequence is: ', Reverse_Complement(seq))
    fasta_file = st.file_uploader("Upload a FASTA or FNA file", type=["fasta", "fna"])
    if fasta_file:
        sequence = parse_fasta(fasta_file)
        if sequence:
            st.write(f'Sequence is: {sequence}')
            # st.write("### Reverse Complement Sequence:")
            st.write(f'Reverse Complement Sequence is: {Reverse_Complement(sequence)}')
    
if st.session_state.page == 'Translation_Table':
    st.title('**Translation Table of a DNA Sequence**')
    st.write('A translation table is a reference used in molecular biology to convert\
              a sequence of nucleotides (codons) into amino acids during the process of translation.\
              This process converts messenger RNA (mRNA) into a protein.')
    seq = st.text_area("Enter Sequence for DNA:", placeholder="Type here...")
    if st.button("Submit"):
        if seq:
            st.write(f'Translation Table of Sequence is: ', Translation_Table(seq)) 

    fasta_file = st.file_uploader("Upload a FASTA or FNA file", type=["fasta", "fna"])
    if fasta_file:
        sequence = parse_fasta(fasta_file)
        if sequence:
            st.write(f'Sequence is: {sequence}')
            # st.write("### Translation Table Sequence:")
            st.write(f'Translation Table of Sequence is: {Translation_Table(sequence)}')

    image = Image.open('Translation_table.png')
    st.image(image, use_column_width=True)

if st.session_state.page == 'Match':
    st.title('**Sequence and Sub-sequence Matcher**')
    st.write('The match function performs a basic search for a sub-sequence within a sequence and\
              returns the starting index of the first match or -1 if no match is found.\
              It compares each substring of seq with sub_seq and returns the index of the first match it encounters.')
    seq = st.text_area("Enter Sequence for DNA:", placeholder="Type here...")
    sub_seq = st.text_input("Enter the sub-sequence to search for:", placeholder="Type here...")
    if st.button("Submit"):
        if seq and sub_seq:
            if Match(seq, sub_seq) != -1:
                st.write(f'Match founded sub-sequence at index:', Match(seq, sub_seq))
            else:
                st.write("No match found.")

    fasta_file = st.file_uploader("Upload a FASTA or FNA file", type=["fasta", "fna"])
    if fasta_file:
        sequence = parse_fasta(fasta_file)
        if sequence:
            st.write(f'Sequence is: {sequence}')
    sub_seq_fasta = st.text_input("Enter the sub-sequence to search in fasta file:", placeholder="Type here...")
    if st.button("Submit Fasta"):
        if fasta_file and sub_seq_fasta:
            if Match(sequence, sub_seq_fasta) != -1:
                st.write(f'Match founded sub-sequence at index:', Match(sequence, sub_seq_fasta))
            else:
                st.write("No match found.")

if st.session_state.page == 'Bad':
    st.title('**Bad Character Heuristic for String Matching**')
    st.write('The Bad Character Heuristic is a technique used in the Boyer-Moore string matching algorithm\
              to optimize the process of searching for a sub-sequence (or pattern) within a main sequence (or text).\
              It helps skip over sections of the main sequence that cannot possibly match the sub-sequence,\
              reducing the number of comparisons needed.')
    seq = st.text_area("Enter Sequence for DNA:", placeholder="Type here...")
    sub_seq = st.text_input("Enter the sub-sequence to search for:", placeholder="Type here...")
    if st.button("Submit"):
        if seq and sub_seq:
            if Badchars(seq, sub_seq) != -1:
                st.write(f'Match founded ({sub_seq}) at index:', Badchars(seq, sub_seq))
            else:
                st.write("No match found.")

    fasta_file = st.file_uploader("Upload a FASTA or FNA file", type=["fasta", "fna"])
    if fasta_file:
        sequence = parse_fasta(fasta_file)
        if sequence:
            st.write(f'Sequence is: {sequence}')
    sub_seq_fasta = st.text_input("Enter the sub-sequence to search in fasta file:", placeholder="Type here...")
    if st.button("Submit Fasta"):
        if fasta_file and sub_seq_fasta:
            if Badchars(sequence, sub_seq_fasta) != -1:
                st.write(f'Match founded sub-sequence at index:', Badchars(sequence, sub_seq_fasta))
            else:
                st.write("No match found.")

if st.session_state.page == 'Index':
    st.title('**Sorted Substrings with Indices**')
    st.write('The IndexSorted function extracts all possible substrings of a given\
              length (ln) from a sequence (seq) and returns them along with their starting positions.\
              It then sorts these substrings lexicographically (alphabetically).')
    seq = st.text_area("Enter Sequence for DNA:", placeholder="Type here...")
    ln = st.number_input("Enter the length of the substring:", min_value=0, max_value=len(seq), step=1)
    ln = int(ln)
    if st.button("Submit"):
        if seq and ln:
            st.write(f'Sorted Substrings with Indices:')
            for substring, index in IndexSorted(seq, ln):
                st.write(f"Substring: '{substring}' at index {index}")

    fasta_file = st.file_uploader("Upload a FASTA or FNA file", type=["fasta", "fna"])
    sequence = ''
    if fasta_file:
        sequence = parse_fasta(fasta_file)
        if sequence:
            st.write(f'Sequence is: {sequence}')
    ln_fasta = st.number_input("Enter the length of the substring to search in fasta file:", min_value=0, max_value=len(sequence), step=1)
    ln_fasta = int(ln_fasta)
    if st.button("Submit Fasta"):
        if fasta_file and ln_fasta:
            st.write(f'Sorted Substrings with Indices:')
            for substring, index in IndexSorted(sequence, ln_fasta):
                st.write(f"Substring: '{substring}' at index {index}")

if st.session_state.page == 'Query':
    st.title("**Pattern Search in Text**")
    st.write("The query function youve shared is designed to perform pattern matching in a string\
              using an index to speed up the search process.\
              It finds the starting positions of the pattern in the string by leveraging a binary search approach.")
    seq = st.text_area("Enter Sequence for DNA:", placeholder="Type here...")
    sub_seq = st.text_input("Enter the sub-sequence to search for:", placeholder="Type here...")
    if st.button("Submit"):
        if seq and sub_seq:
            result = query(seq, sub_seq, IndexSorted(seq, len(sub_seq)))
            if result:
                st.write(f"Pattern found at positions: {result}")
            else:
                st.write("Pattern not found.")

    fasta_file = st.file_uploader("Upload a FASTA or FNA file", type=["fasta", "fna"])
    if fasta_file:
        sequence = parse_fasta(fasta_file)
        if sequence:
            st.write(f'Sequence is: {sequence}')
    sub_seq_fasta = st.text_input("Enter the sub-sequence to search in fasta file:", placeholder="Type here...")
    if st.button("Submit Fasta"):
        if fasta_file and sub_seq_fasta:
            result = query(sequence, sub_seq_fasta, IndexSorted(sequence, len(sub_seq_fasta)))
            if result:
                st.write(f"Pattern found at positions: {result}")
            else:
                st.write("Pattern not found.")

if st.session_state.page == 'Suffix_array_construction':
    st.title("**Suffix array construction**")
    st.write('The general rule is that, in iteration i, all suffixes are sorted \
              according to their first 2ⁱ characters only.')
    st.write('$ : 0, A : 1, C : 2, G : 3, T : 4.')
    seq = st.text_area("Enter Sequence for DNA:", placeholder="Type here...")
    if st.button("Submit"):
        if seq:
            st.write('Suffix array:')
            iterations = Suffix_array_construction(seq)
            for iteration in iterations:
                st.write(iteration)
    
    fasta_file = st.file_uploader("Upload a FASTA or FNA file", type=["fasta", "fna"])
    if fasta_file:
        sequence = parse_fasta(fasta_file)
        if sequence:
            st.write(f'Sequence is: {sequence}')
            st.write('Suffix array:')
            iterations_fasta = Suffix_array_construction(sequence)
            for iteration_fasta in iterations_fasta:
                st.write(iteration_fasta)

if st.session_state.page == 'Overlap':
    st.title('**Overlap**')
    st.write('The overlap function is used to find the length of the overlapping part between two strings, \
             a and b, where the overlap is defined as the end of string a matching the beginning of string b.')
    seq1 = st.text_area("Enter Sequence for DNA:", placeholder="Type here...")
    seq2 = st.text_area("Enter Sequence2 for DNA to Overlap with Sequence1:", placeholder="Type here...")
    if st.button("Submit"):
        if seq1 and seq2:
            if overlap(seq1, seq2) != 0:
                st.write(f'Overlap with K-:',overlap(seq1, seq2))
            else:
                st.write("No Overlap.")

if st.session_state.page == 'native_overlap':
    st.title('**Native Overlap**')
    st.write('The Native Overlap algorithm determines the largest overlap between pairs of DNA sequences in a list')
    seq = st.text_area("Enter Sequence for DNA:", placeholder="Type here...")
    k = st.number_input("Enter the minimum length of the overlap:", min_value=0, step=1)
    k = int(k)
    if st.button("Submit"):
        if seq and k:
            result = native_overlap(seq, k)
            st.subheader("**Overlaps Found:**")
            for key, value in result.items():
                st.write(f"{key} : {value}")
    

