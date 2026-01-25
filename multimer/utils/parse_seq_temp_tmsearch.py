#!/usr/bin/env python3
"""
Template CSV/TSV Alignment Parser

Extracts template alignments from CSV/TSV files with alignment data and saves to standardized CSV format.
Parses aligned sequences to get precise residue-to-residue mappings.
"""

import csv
import pandas as pd
import argparse
from pathlib import Path
from typing import List, Tuple, Dict, Optional


class TemplateCSVParser:
    def __init__(self):
        self.templates = []
    
    def parse_template_file(self, template_file: str, delimiter: str = None) -> List[Dict]:
        """Parse template alignment file (CSV/TSV) and extract alignments."""
        
        # Auto-detect delimiter if not specified
        if delimiter is None:
            with open(template_file, 'r') as f:
                first_line = f.readline()
                if '\t' in first_line:
                    delimiter = '\t'
                else:
                    delimiter = ','
        
        # Read the file
        try:
            df = pd.read_csv(template_file, delimiter=delimiter)
        except Exception as e:
            print(f"Error reading file: {e}")
            return []
        
        
        templates = []
        for idx, row in df.iterrows():
            template_data = self._parse_alignment_row(row)
            if template_data:
                templates.append(template_data)
        
        return templates
    
    def _parse_alignment_row(self, row: pd.Series) -> Optional[Dict]:
        """Parse individual alignment row to extract alignment data."""
        try:
            # Extract basic info
            template_name = str(row['target'])
            query_seq = str(row['qaln'])
            template_seq = str(row['taln'])
            query_start = int(row['qstart'])
            template_start = int(row['tstart'])
            
            # Extract aligned positions
            query_indices, template_indices = self._extract_aligned_positions(
                query_seq, template_seq, query_start, template_start
            )
            
            if not query_indices or not template_indices:
                return None
            
            return {
                'template_name': template_name.replace("_",""),
                'query_indices': query_indices,
                'template_indices': template_indices,
                'alignment_length': len(query_indices),
                'query_start': query_start,
                'template_start': template_start,
                'sum_probs': row.get('sum_probs', None)
            }
            
        except (KeyError, ValueError) as e:
            print(f"Error parsing row: {e}")
            return None
    
    def _extract_aligned_positions(self, query_seq: str, template_seq: str, 
                                 query_start: int, template_start: int) -> Tuple[List[int], List[int]]:
        """Extract aligned residue positions from sequence alignment."""
        query_indices = []
        template_indices = []
        
        query_pos = query_start
        template_pos = template_start
        
        for q_char, t_char in zip(query_seq, template_seq):
            # print()
            # Both positions have residues (aligned)
            
            if not q_char == "-" and not t_char == "-":
                query_indices.append(query_pos)
                template_indices.append(template_pos)
                template_pos += 1
                query_pos += 1
        
        return query_indices, template_indices
    
    def save_to_csv(self, templates: List[Dict], output_file: str):
        """Save extracted template data to CSV file."""
        with open(output_file, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            
            # # Write header
            # writer.writerow(['template_name', 'query_indices', 'template_indices', 
            #                'alignment_length', 'query_start', 'template_start', 'sum_probs'])

            writer.writerow(['template_name', 'query_indices', 'template_indices'])
            
            # Write data
            for template in templates:
                query_indices_str = ','.join(map(str, template['query_indices']))
                template_indices_str = ','.join(map(str, template['template_indices']))
                query_indices = [i for i in template['query_indices']]
                template_indices = [i for i in template['template_indices']]
                
                writer.writerow([
                    template['template_name'],
                    query_indices,
                    template_indices
                    # template['alignment_length'],
                    # template['query_start'],
                    # template['template_start'],
                    # template.get('sum_probs', '')
                ])
    
    def save_to_detailed_csv(self, templates: List[Dict], output_file: str):
        """Save to CSV with one alignment pair per row (more detailed format)."""
        with open(output_file, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            
            # Write header
            writer.writerow(['template_name', 'query_index', 'template_index', 
                           'alignment_position', 'sum_probs'])
            
            # Write data
            for template in templates:
                template_name = template['template_name']
                sum_probs = template.get('sum_probs', '')
                
                for pos, (q_idx, t_idx) in enumerate(zip(template['query_indices'], 
                                                        template['template_indices'])):
                    writer.writerow([template_name, q_idx, t_idx, pos + 1, sum_probs])
    
    def save_to_alphafold_format(self, templates: List[Dict], output_file: str):
        """Save in a format specifically useful for AlphaFold3 template processing."""
        with open(output_file, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            
            # # Write header with AlphaFold3-relevant fields
            # writer.writerow(['template_name', 'query_indices', 'template_indices', 
            #                'num_aligned_residues', 'query_range', 'template_range', 
            #                'alignment_score'])
            writer.writerow(['template_name', 'query_indices', 'template_indices'])
            
            # Write data
            for template in templates:
                query_indices = template['query_indices']
                template_indices = template['template_indices']
                
                query_indices_str = ','.join(map(str, query_indices))
                template_indices_str = ','.join(map(str, template_indices))
                
                query_range = f"{min(query_indices)}-{max(query_indices)}"
                template_range = f"{min(template_indices)}-{max(template_indices)}"
                
                writer.writerow([
                    template['template_name'],
                    query_indices_str,
                    template_indices_str,
                    # len(query_indices),
                    # query_range,
                    # template_range,
                    # template.get('sum_probs', '')
                ])
    
    def print_summary(self, templates: List[Dict]):
        """Print summary statistics of parsed templates."""
        if not templates:
            print("No templates found.")
            return
        
        print(f"\nTemplate Summary:")
        print(f"Total templates: {len(templates)}")
        
        # Alignment length statistics
        lengths = [t['alignment_length'] for t in templates]
        print(f"Alignment lengths: min={min(lengths)}, max={max(lengths)}, avg={sum(lengths)/len(lengths):.1f}")
        
        # Template names
        template_names = [t['template_name'] for t in templates]
        unique_templates = set(template_names)
        print(f"Unique templates: {len(unique_templates)}")
        
        # Score statistics if available
        scores = [t.get('sum_probs') for t in templates if t.get('sum_probs') is not None]
        if scores:
            print(f"Alignment scores: min={min(scores):.1f}, max={max(scores):.1f}, avg={sum(scores)/len(scores):.1f}")
        
        print(f"\nTop 5 templates by alignment length:")
        sorted_templates = sorted(templates, key=lambda x: x['alignment_length'], reverse=True)
        for i, template in enumerate(sorted_templates[:5]):
            print(f"  {i+1}. {template['template_name']}: {template['alignment_length']} residues")

def run_seq_temp_parser(
    template_file,
    output='template_alignments.csv',
    delimiter='auto',
    detailed=False,
    alphafold=False,
    summary=False,
    verbose=False
):
    """
    Parse a template alignment file (CSV/TSV) and save it in the requested format.

    Args:
        template_file (str): Input file path
        output (str): Output CSV file path
        delimiter (str): ',' | '\t' | 'auto'
        detailed (bool): Save one pair per row
        alphafold (bool): Save in AlphaFold3-optimized format
        summary (bool): Print summary statistics
        verbose (bool): Verbose output
    Returns:
        int: 0 on success, 1 on error
    """
    if not Path(template_file).exists():
        print(f"Error: Template file '{template_file}' not found.")
        return 1

    delim = None if delimiter == 'auto' else delimiter

    # from your_module import TemplateCSVParser  # adjust import if needed
    parser_obj = TemplateCSVParser()

    if verbose:
        print(f"Parsing template file: {template_file}")

    templates = parser_obj.parse_template_file(template_file, delim)

    if not templates:
        print("Warning: No template alignments found in file.")
        return 1

    if verbose or summary:
        parser_obj.print_summary(templates)

    if alphafold:
        parser_obj.save_to_alphafold_format(templates, output)
        # print(f"AlphaFold3-format template alignments saved to: {output}")
    elif detailed:
        parser_obj.save_to_detailed_csv(templates, output)
        # print(f"Detailed template alignments saved to: {output}")
    else:
        parser_obj.save_to_csv(templates, output)
        # print(f"Template alignments saved to: {output}")

    return 0
