#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Author: Tomasz Stawski, tomasz.stawski@bam.de
"""
Script to apply author headers, British spelling, and em-dash removal across all Python scripts.
Neutral language without excessive adjectives.
"""

import os
import glob
import re

AUTHOR_LINE = "# Author: Tomasz Stawski, tomasz.stawski@bam.de\n"
AUTHOR_STR = "Tomasz Stawski, tomasz.stawski@bam.de"

BRITISH_REPLACEMENTS = [
    (r'\bmodeling\b', 'modelling'),
    (r'\bModeling\b', 'Modelling'),
    (r'\bmodeled\b', 'modelled'),
    (r'\bModeled\b', 'Modelled'),
    (r'\bbehavior\b', 'behaviour'),
    (r'\bBehavior\b', 'Behaviour'),
    (r'\bbehaviors\b', 'behaviours'),
    (r'\bBehaviors\b', 'Behaviours'),
    (r'\banalyze\b', 'analyse'),
    (r'\bAnalyze\b', 'Analyse'),
    (r'\banalyzed\b', 'analysed'),
    (r'\bAnalyzed\b', 'Analysed'),
    (r'\banalyzing\b', 'analysing'),
    (r'\bAnalyzing\b', 'Analysing'),
    (r'\bcharacterize\b', 'characterise'),
    (r'\bCharacterize\b', 'Characterise'),
    (r'\bcharacterized\b', 'characterised'),
    (r'\bCharacterized\b', 'Characterised'),
    (r'\bcharacterizing\b', 'characterising'),
    (r'\bCharacterizing\b', 'Characterising'),
    (r'\bcharacterization\b', 'characterisation'),
    (r'\bCharacterization\b', 'Characterisation'),
    (r'\bsynthesize\b', 'synthesise'),
    (r'\bSynthesize\b', 'Synthesise'),
    (r'\bsynthesized\b', 'synthesised'),
    (r'\bSynthesized\b', 'Synthesised'),
    (r'\bemphasize\b', 'emphasise'),
    (r'\bEmphasize\b', 'Emphasise'),
    (r'\bemphasized\b', 'emphasised'),
    (r'\bEmphasized\b', 'Emphasised'),
]

def process_file(filepath):
    with open(filepath, 'r', encoding='utf-8') as f:
        content = f.read()

    modified = False
    
    # 1. Replace em-dash characters
    if ':' in content:
        # Check context
        content = content.replace("N/A (Induction)", "N/A (Induction)")
        content = content.replace(":", ":")
        modified = True

    # 2. British spelling in comments, docstrings, and strings
    # We do safe line-by-line replacement for lines that are comments, docstrings, or string literals,
    # avoiding keywords in code like 'scipy.optimize' or 'from lmfit import minimize' or 'color='
    lines = content.split('\n')
    new_lines = []
    
    in_docstring = False
    docstring_char = None
    
    for line in lines:
        stripped = line.strip()
        
        # Check docstring toggles
        if '"""' in line:
            # count occurrences
            count = line.count('"""')
            if count % 2 == 1:
                in_docstring = not in_docstring
        elif "'''" in line:
            count = line.count("'''")
            if count % 2 == 1:
                in_docstring = not in_docstring
                
        is_comment = stripped.startswith('#')
        is_string_or_doc = in_docstring or is_comment or ('"' in line or "'" in line)
        
        if is_string_or_doc:
            # Check replacements
            line_orig = line
            for pat, repl in BRITISH_REPLACEMENTS:
                # Do not replace if it's an import or code definition
                if 'import' in line and not is_comment and not in_docstring:
                    continue
                line = re.sub(pat, repl, line)
            if line != line_orig:
                modified = True
                
        new_lines.append(line)
        
    content = '\n'.join(new_lines)

    # 3. Add author if not present
    if AUTHOR_STR not in content:
        lines = content.split('\n')
        insert_idx = 0
        if len(lines) > 0 and lines[0].startswith('#!'):
            insert_idx = 1
        if len(lines) > insert_idx and 'coding:' in lines[insert_idx]:
            insert_idx += 1
            
        lines.insert(insert_idx, AUTHOR_LINE.strip())
        content = '\n'.join(lines)
        modified = True

    if modified:
        with open(filepath, 'w', encoding='utf-8') as f:
            f.write(content)
        return True
    return False


def main():
    root_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
    py_files = sorted(glob.glob(os.path.join(root_dir, '**', '*.py'), recursive=True))
    
    updated = 0
    for pf in py_files:
        if process_file(pf):
            updated += 1
            print(f"Updated: {os.path.relpath(pf, root_dir)}")

    print(f"\nTotal Python files processed: {len(py_files)}")
    print(f"Total Python files updated: {updated}")

if __name__ == '__main__':
    main()
