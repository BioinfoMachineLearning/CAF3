#!/usr/bin/env python3
"""
HHR Template Alignment Parser

Extracts template alignments from HHR files and saves to CSV format.
Parses the full alignment blocks to get precise residue-to-residue mappings.
"""

import re
import csv
import argparse
from pathlib import Path
from typing import List, Tuple, Dict, Optional


class HHRParser:
    def __init__(self):
        self.templates = []
    
    def parse_hhr_file(self, hhr_file: str) -> List[Dict]:
        """Parse HHR file and extract template alignments."""
        with open(hhr_file, 'r') as f:
            content = f.read()
        
        # Debug: Print first few lines to understand structure
        lines = content.split('\n')
        template_blocks = re.split(r'\n>', content)
        
        if len(template_blocks) <= 1:
            template_blocks = re.split(r'\nNo\s+\d+', content)
        
        
        templates = []
        for i, block in enumerate(template_blocks):
            if i == 0 and not block.strip().startswith('>'):
                # First block might be header
                continue
            
            template_data = self._parse_template_block(block)
            if template_data:
                templates.append(template_data)
            else:
                print(f"DEBUG: Failed to parse block {i}")
        
        return templates
    
    def _parse_template_block(self, block: str) -> Optional[Dict]:
        """Parse individual template block to extract alignment data."""
        lines = block.strip().split('\n')
        if not lines:
            return None
        
        # Extract template name from first line
        template_name = lines[0].split()[0]
        
        # Find alignment section - look for the actual alignment pattern
        alignment_start = None
        for i, line in enumerate(lines):
            # Look for lines that start alignment blocks (usually have "No " followed by alignment)
            # or look for query/template alignment patterns
            if (line.strip().startswith('No ') and 
                i + 1 < len(lines) and 
                ('Q ' in lines[i + 1] or 'T ' in lines[i + 1])):
                alignment_start = i + 1
                break
            # Alternative: look directly for Q/T alignment patterns
            elif line.strip().startswith('Q ') and len(line.split()) >= 4:
                alignment_start = i
                break
        
        if alignment_start is None:
            # Try a more flexible approach - look for any line with alignment pattern
            for i, line in enumerate(lines):
                if self._is_alignment_line(line):
                    alignment_start = i
                    break
        
        if alignment_start is None:
            return None
        
        # Parse alignment blocks
        query_indices, template_indices = self._parse_alignment_blocks(
            lines[alignment_start:]
        )
        
        if not query_indices or not template_indices:
            return None
        
        return {
            'template_name': template_name.replace("_",""),
            'query_indices': query_indices,
            'template_indices': template_indices
        }
    
    def _is_alignment_line(self, line: str) -> bool:
        """Check if line looks like an alignment line."""
        parts = line.strip().split()
        if len(parts) < 4:
            return False
        
        # Check for Q/T alignment pattern: "Q identifier start_pos sequence end_pos"
        if parts[0] in ['Q', 'T']:
            try:
                int(parts[2])  # start position should be integer
                return True
            except ValueError:
                pass
        
        return False
    
    def _parse_alignment_blocks(self, alignment_lines: List[str]) -> Tuple[List[int], List[int]]:
        """Parse alignment blocks to extract aligned residue indices."""
        query_indices = []
        template_indices = []
        
        i = 0
        while i < len(alignment_lines):
            line = alignment_lines[i].strip()
            
            # Look for query line (starts with 'Q ')
            if line.startswith('Q ') and len(line.split()) >= 4:
                query_parts = line.split()
                try:
                    query_start = int(query_parts[2])
                    query_seq = query_parts[3]
                    query_end = int(query_parts[4]) if len(query_parts) >= 5 else None
                except (ValueError, IndexError):
                    i += 1
                    continue
                
                # Look for corresponding template line (starts with 'T ')
                template_line = None
                template_idx = None
                for j in range(i + 1, min(i + 5, len(alignment_lines))):  # Look within next 4 lines
                    if alignment_lines[j].strip().startswith('T ') and "Consensus" not in alignment_lines[j].strip() and "ss_pred" not in alignment_lines[j].strip():
                        template_line = alignment_lines[j].strip()
                        template_idx = j
                        break
                
                if template_line:
                    template_parts = template_line.split()
                    if len(template_parts) >= 4 :
                        try:
                            template_start = int(template_parts[2])
                            template_seq = template_parts[3]
                            
                            # Extract aligned positions
                            q_indices, t_indices = self._extract_aligned_positions(
                                query_seq, template_seq, query_start, template_start
                            )
                            query_indices.extend(q_indices)
                            template_indices.extend(t_indices)
                            
                            # Skip to after the template line
                            i = template_idx + 1
                            continue
                        except (ValueError, IndexError):
                            pass
            
            i += 1
        
        return query_indices, template_indices
    
    def _extract_aligned_positions(self, query_seq: str, template_seq: str, 
                                 query_start: int, template_start: int) -> Tuple[List[int], List[int]]:
        """Extract aligned residue positions from sequence alignment."""
        query_indices = []
        template_indices = []
        
        query_pos = query_start
        template_pos = template_start
        
        for q_char, t_char in zip(query_seq, template_seq):

            if not q_char == "-" and not t_char == "-":
                query_indices.append(query_pos)
                template_indices.append(template_pos)
            if q_char != "-":
                query_pos += 1
            if t_char != "-":
                template_pos += 1

        
        return query_indices, template_indices
    
    def save_to_csv(self, templates: List[Dict], output_file: str):
        """Save extracted template data to CSV file."""
        with open(output_file, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            
            # Write header
            writer.writerow(['template_name', 'query_indices', 'template_indices'])
            
            # Write data
            for template in templates:
                query_indices_str = ','.join(map(str, template['query_indices']))
                template_indices_str = ','.join(map(str, template['template_indices']))
                query_indices = [i - 1 for i in template['query_indices']]
                template_indices = [i - 1 for i in template['template_indices']]
                
                writer.writerow([
                    template['template_name'],
                    query_indices,
                    template_indices
                ])
    
    def save_to_detailed_csv(self, templates: List[Dict], output_file: str):
        """Save to CSV with one alignment pair per row (more detailed format)."""
        with open(output_file, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            
            # Write header
            writer.writerow(['template_name', 'query_index', 'template_index'])
            
            # Write data
            for template in templates:
                template_name = template['template_name']
                for q_idx, t_idx in zip(template['query_indices'], template['template_indices']):
                    writer.writerow([template_name, q_idx, t_idx])

def run_hhr_parser(
    hhr_file,
    output='template_alignments.csv',
    detailed=False,
    verbose=False
):
    """
    Parse an HHR file and save template alignments to CSV.

    Args:
        hhr_file (str): Path to input .hhr file
        output (str): Path to output CSV file
        detailed (bool): If True, save one aligned pair per row
        verbose (bool): If True, print progress
    Returns:
        int: 0 on success, 1 on failure
    """
    # Check if input file exists
    if not Path(hhr_file).exists():
        print(f"Error: HHR file '{hhr_file}' not found.")
        return 1

    # Import your parser class here
    # from your_module import HHRParser  # change this if needed
    hhr_parser = HHRParser()

    if verbose:
        print(f"Parsing HHR file: {hhr_file}")

    templates = hhr_parser.parse_hhr_file(hhr_file)

    if not templates:
        print("Warning: No template alignments found in HHR file.")
        return 1

    # if verbose:
    #     # print(f"Found {len(templates)} template alignments")
    #     for template in templates:
    #         print(f"  {template['template_name']}: {len(template['query_indices'])} aligned residues")

    # Save to CSV
    if detailed:
        hhr_parser.save_to_detailed_csv(templates, output)
    else:
        hhr_parser.save_to_csv(templates, output)

    # print(f"Template alignments saved to: {output}")
    return 0
