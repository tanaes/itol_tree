#!/usr/bin/env python

import warnings
import palettable
import requests
import os
import numpy as np
import pandas as pd
from zipfile import ZipFile
from ete3 import Tree
from random import randint

### iToL metadata datafile routines

def get_qualitative_colors(fields, palette=palettable.colorbrewer.qualitative.Set3_12):
    try:
        colors = [palette.hex_colors[x] for x in range(len(fields))]
    except IndexError:
        warnings.warn("Too many fields for palette. Recycling colors")
        colors = [palette.hex_colors[x % palette.number] for x in range(len(fields))]
    return(colors)


def format_iTOL_multibar(fields, metadata, 
                         tree_ref_col=None, 
                         field_colors=None, 
                         field_labels=None, 
                         dataset_label='Multibar Chart', 
                         dataset_color=None, 
                         legend=True, 
                         unstacked=False,
                         width=1000, 
                         margin=0
                        ):
    """
    fields: array of columns titles in [metadata] to chart
    metadata: pd.df containing samples to graph and data
    """
    
    if field_labels is None:
        field_labels=fields
    
    if field_colors is None:
        field_colors=get_qualitative_colors(fields)
    
    if tree_ref_col is None:
        tree_ref_col=metadata.columns[0]

    if dataset_color is None:
        dataset_color = "#%06x" % randint(0, 0xFFFFFF)
    
    outstring = ''
    
    outstring += 'DATASET_MULTIBAR\n'
    outstring += 'SEPARATOR TAB\n'
    outstring += 'DATASET_LABEL\t%s\n' % dataset_label
    outstring += 'COLOR\t%s\n' % dataset_color
    outstring += 'FIELD_COLORS\t%s\n' % '\t'.join(field_colors)
    outstring += 'FIELD_LABELS\t%s\n' % '\t'.join(field_labels)    
    
    if legend:
        outstring += 'LEGEND_TITLE\tDataset legend\n'
        outstring += 'LEGEND_SHAPES\t%s\n' % '\t'.join(['1']*len(fields))
        outstring += 'LEGEND_COLORS\t%s\n' % '\t'.join(field_colors)
        outstring += 'LEGEND_LABELS\t%s\n' % '\t'.join(field_labels)
    
    outstring += 'MARGIN\t%s\n' % margin
    outstring += 'WIDTH\t%s\n' % width

    if unstacked:
        outstring += 'ALIGN_FIELDS\t1\n'
        
    outstring += 'DATA\n'
    for index, row in metadata.iterrows():
        outstring += row[tree_ref_col] + '\t%s\n' % '\t'.join([str(row[x]) for x in fields])

    return(outstring)


def format_iTOL_heatmap(fields, metadata, 
                         field_colors=None, 
                         field_labels=None,
                         field_tree=None,
                         show_field_tree=False,
                         color_minmax = ['#ff0000','#0000ff'],
                         mid_color = '#ffffff',
                         color_limits = False,
                         dataset_label='Heatmap', 
                         dataset_color=None, 
                         legend=False, 
                         strip_width=2,
                         margin=0
                        ):
    """
    fields: array of columns titles in [metadata] to chart
    metadata: pd.df containing samples to graph and data
    """
    
    if field_labels is None:
        field_labels=fields
    
    if dataset_color is None:
        dataset_color = "#%06x" % randint(0, 0xFFFFFF)

    outstring = ''
    
    outstring += 'DATASET_HEATMAP\n'
    outstring += 'SEPARATOR TAB\n'
    outstring += 'DATASET_LABEL\t%s\n' % dataset_label
    outstring += 'COLOR\t%s\n' % dataset_color
    outstring += 'FIELD_LABELS\t%s\n' % '\t'.join(field_labels)    
    
    if field_tree is not None:
        outstring += 'FIELD_TREE\t%s\n' % field_tree.write()
        
    if legend:
        outstring += 'LEGEND_TITLE\tDataset legend\n'
        outstring += 'LEGEND_SHAPES\t%s\n' % '\t'.join(['1']*len(fields))
        outstring += 'LEGEND_COLORS\t%s\n' % '\t'.join(field_colors)
        outstring += 'LEGEND_LABELS\t%s\n' % '\t'.join(field_labels)

    outstring += 'MARGIN\t%s\n' % margin

    if strip_width is not None:
        outstring += 'STRIP_WIDTH\t%s\n' % strip_width
    
    outstring += 'COLOR_MIN\t%s\n' % color_minmax[0]
    outstring += 'COLOR_MAX\t%s\n' % color_minmax[1]
    
    if mid_color:
        outstring += 'USE_MID_COLOR\t1\n'
        outstring += 'COLOR_MID\t%s\n' % mid_color
    
    if color_limits:
            # By default, color gradients will be calculated based on dataset values.
            # You can force different values to use in the calculation by setting the values below:
            # provide as three int list, min mid max

            outstring += 'USER_MIN_VALUE\t%s\n' % color_limits[0]
            outstring += 'USER_MID_VALUE\t%s\n' % color_limits[1]
            outstring += 'USER_MAX_VALUE\t%s\n' % color_limits[2]

    outstring += 'DATA\n'
    for index, row in metadata.iterrows():
        outstring += index + '\t%s\n' % '\t'.join([str(row[x]) for x in fields])

    return(outstring)


def get_defining_tips(t):
    first = True
    for node in t.traverse("postorder"):
        if node.is_leaf():
            if first:
                node1 = node.name
                first = False
            else:
                node2 = node.name
        # Do some analysis on node
        # print node.name
    
    return(node1, node2)


def format_itol_range_colors(tree, metadata, color_cat, tree_ref_col=None, color_vals=None, cat_colors=None, cat_labels=None):
    
    if color_vals is None:
        color_vals = list(set(metadata[color_cat]))
        
    if cat_colors is None:
        cat_colors = get_qualitative_colors(color_vals, palette=palettable.colorbrewer.qualitative.Set3_12)
    
    if cat_labels is None:
        cat_labels = color_vals
    
    if tree_ref_col is None:
        tree_ref_col=metadata.columns[0]
        
    # pick out values in color_cat
    outstring = 'TREE_COLORS\n'
    outstring += 'SEPARATOR TAB\n'
    outstring += 'DATA\n'    
    
    # for each value, get the defining pair of tips
    for val in color_vals:
        i = color_vals.index(val)
        
        defining_tips = get_defining_tips(
                         tree.get_common_ancestor(
                          metadata.loc[metadata[color_cat] == val,
                                       tree_ref_col].values.tolist()))
        
        outstring += '{0}\t{1}\t{2}\t{3}\n'.format('|'.join(defining_tips),
                                                   'range',
                                                   cat_colors[i],
                                                   cat_labels[i])
        
    return(outstring)


def format_itol_tip_colors(metadata, color_cat, tree_ref_col=None, color_vals=None, cat_colors=None, cat_labels=None):
    
    if color_vals is None:
        color_vals = list(set(metadata[color_cat]))
        
    if cat_colors is None:
        cat_colors = get_qualitative_colors(color_vals, palette=palettable.colorbrewer.qualitative.Set3_12)
    
    if cat_labels is None:
        cat_labels = color_vals
    
    if tree_ref_col is None:
        tree_ref_col=metadata.columns[0]
        
    # pick out values in color_cat
    outstring = 'TREE_COLORS\n'
    outstring += 'SEPARATOR TAB\n'
    outstring += 'DATA\n'    
    
    # for each value, get the defining pair of tips
    for val in color_vals:
        i = color_vals.index(val)
        
        tips = metadata.loc[metadata[color_cat] == val,
                                       tree_ref_col].values.tolist()
        # #8015,label,#0000ff
        for tip in tips:
            outstring += '{0}\t{1}\t{2}\n'.format(tip,
                                                  'label',
                                                   cat_colors[i])
        
    return(outstring)

def format_iTOL_simplebar(field, metadata, 
                         tree_ref_col=None, 
                         bar_color=None, 
                         dataset_label='Simplebar Chart', 
                         dataset_color="#%06x" % randint(0, 0xFFFFFF),
                         width=1000,
                         margin=0
                        ):
    """
    fields: array of columns titles in [metadata] to chart
    metadata: pd.df containing samples to graph and data
    """
    
    if bar_color is None:
        bar_color = '#ff0000'
    
    if tree_ref_col is None:
        tree_ref_col=metadata.columns[0]
    
    outstring = ''
    
    outstring += 'DATASET_SIMPLEBAR\n'
    outstring += 'SEPARATOR TAB\n'
    outstring += 'DATASET_LABEL\t%s\n' % dataset_label
    outstring += 'COLOR\t%s\n' % dataset_color    

    outstring += 'MARGIN\t%s\n' % margin

    outstring += 'WIDTH\t%s\n' % width
        
    outstring += 'DATA\n'
    for index, row in metadata.iterrows():
        outstring += row[tree_ref_col] + '\t%s\n' % row[field]

    return(outstring)


### iToL batch access API routines

def prepare_itol_files(tree, dataset_dict = {}, outdir='./', basename='iToL_tree'):

    tree_fp = os.path.join(outdir, '%s.tree' % basename)

    tree.write(format=1, outfile=tree_fp)

    zip_files = []

    zip_files += [tree_fp]
    
    for i, dataset in enumerate(sorted(dataset_dict)):
        data_fp = os.path.join(outdir, '{0}_{1}_{2}.txt'.format(basename, dataset, i))
        with open(data_fp, 'w') as f:
            f.write(dataset_dict[dataset])

        zip_files += [data_fp]

    zip_fp = os.path.join(outdir, '%s.zip' % basename)

    with ZipFile(zip_fp, 'w') as datazip:
        for fp in zip_files:
            datazip.write(fp)
    
    return(zip_fp)


# Jon's Upload ID: fRTsUmgrGbFq1jNuIwaZhQ

def upload_itol_zip(zip_fp, uploadID, projectName, treeName, itol_upload_url = "http://itol.embl.de/batch_uploader.cgi"):

    files = {'zipFile': open(zip_fp,'rb')}
    values = {'uploadID': uploadID, 'projectName': projectName, 'treeName': treeName}

    r_u = requests.post(itol_upload_url, files=files, data=values)

    if r_u.text.startswith('SUCCESS'):
        tree_num = r_u.text.strip().split(' ')[1]
        return(r_u, tree_num)
    else:
        print('Problem uploading tree! \n %s' % r_u.text)
        return(r_u, None)


def download_itol_tree(tree_num, img_fp, file_format='pdf', optional_params={}, itol_download_url = 'http://itol.embl.de/batch_downloader.cgi'):
    # dl_values = {
    # 'tree': str(tree_num),
    # 'format': 'pdf',
    # 'range_mode': '1', # Colored ranges display style. Possible values: 0,1 or 2 (0=off, 1=cover labels only, 2=cover full clades)<
    # 'include_ranges_legend': '1', # Include colored ranges legend. Possible values: 0 or 1
    # 'display_mode': '2', # (1=normal, 2=circular, 3=unrooted)
    # 'leaf_sorting': '1', # Possible values: 1 or 2 (1=normal sorting, 2=no sorting)
    # 'label_display': '1', # Possible values: 0 or 1 (0=hide labels, 1=show labels)
    # 'align_labels': '1', # Possible values: 0 or 1 (0=labels not aligned, 1=labels aligned)
    # 'ignore_branch_length': '0', # Possible values: 0 or 1
    # 'datasets_visible': '0,1,2' # Comma delimited list of datasets to display (starting with 0, e.g. datasets_visible=0,2,5)
    # }

    params = optional_params

    params['tree'] = tree_num
    params['format'] = file_format

    r_d = requests.post(itol_download_url, data = params)

    with open(img_fp, "wb") as tree_f:
        tree_f.write(r_d.content)

    return(r_d)



