#!/usr/bin/env python3

'''
POSSIBLE PROBLEMS

- If you don't use the FastRelax cycle to improve the sequences, the names of the files change. In this case, you should remove the cycle_1 from the csv of Rosetta and CUTRE
This could have an easy fix

'''

'''
things to do:

- Creating files similar to johnbercow for synthesis
- Implementing johnbercow for our use in this webapp
- Pointing out hotspot residue in the PDB viewer from molstar
- Creating campaigns from the webapp
'''

import os
import pandas as pd
import dash
import dash_molstar as dashmolstar
from dash_molstar.utils import molstar_helper
from dash import dash_table
import dash_bio as dashbio
from dash_molstar.utils.representations import Representation
from dash import html
import dash_bootstrap_components as dbc
from dash import Dash, dcc, Input, Output, State, callback
from dash.exceptions import PreventUpdate
import dash_bio.utils.ngl_parser as ngl_parser
import plotly.express as px
import subprocess
import re
from Bio import PDB
import argparse
import socket
import numpy as np
import time
import os.path
from BFmonitor.utils.hits_utils import *
from BFmonitor.utils.generic_utils import *
from BFmonitor.utils.plotting_utils import *
from BFmonitor.utils.input_utils import *
import base64
import io
import tempfile
from dash_molstar.utils.target import Target



# Parsing port number and host
parser = argparse.ArgumentParser()
parser.add_argument("--port", "-p", help = "choose port to run the webapp")
parser.add_argument("--host", "-host", help = "choose host to run the webapp")
parser.add_argument("--debug", "-d", help = "launch app in debug mode")

args, unknown = parser.parse_known_args()

# Set localhost and 8051 as host and port by default
if not args.port: port_number = 8050
else: port_number = args.port
if not args.host: hostname = socket.gethostname()
else: hostname = args.host 
if not args.debug: debug_mode = False
else: debug_mode = True

# Messsssssinesssss
working_dir = os.getcwd()
output_dir = os.path.join(working_dir, 'output')
directories_list = get_working_directories(working_dir)
if not directories_list:
    directories_list=[working_dir]
designs_list=[]

script_dir = os.path.dirname(os.path.abspath(__file__))
#Initial lists for extraction dropdowns

initial_organisms=[
    "Arabidopsis thaliana",
    "Bacillus subtilis",
    "Caenorhabditis elegans",
    "Chlamydomonas reinhardtii",
    "Danio rerio",
    "Drosophila melanogaster",
    "Homo sapiens",
    "Mus musculus",
    "Nicotiana tabacum",
    "Pseudomonas putida",
    "Saccharomyces cerevisiae",
    "Escherichia coli general",
]

# Styles
bt_style = {"align-items": "center", "background-color": "#F2F3F4", "border": "2px solid #000",
            "box-sizing": "border-box", "color": "#000", "cursor": "pointer", "display": "inline-flex",
             'padding':'0.3em 1.2em', 'margin':'0.5em 0.3em 0.3em 0',
            "font-size": "0.9em", 'font-weight':'500', 'border-radius':'2em', 'text-align':'center',
            'transition':'all 0.2s'}
dropdown_style={'background-color':'#F2F3F4',
                "color": "#000", "cursor": "pointer",'width':'260px',
                 'margin':'0.5em 0.3em 0.3em 0',
                "font-size": "1.2em", 'font-weight':'500', 'text-align':'center',
                'transition':'all 0.2s'}

table_viewer_div_style={'display': 'flex', 'justifyContent': 'flex-start'}
table_div_style={'marginRight': '30px'}
table_style={'overflowX': 'auto','width': '100%', 'margin': '0 auto'}
table_cell_style={'minWidth': '150px', 'width': '150px', 'maxWidth': '150px', 'overflow': 'hidden', 'textOverflow': 'ellipsis', 'padding': '5px',}
title_style = {"margin-left": "15px", "margin-top": "15px", "margin-bottom": "0em", "color": "Black", "font-family" : "Helvetica", "font-size":"2.5em"}
box_style3 = {"font-size":"0.9em",'padding':'0.3em 1.2em','width':"18%", "margin-left": "0%","margin-right": "1%", "color": "black","font-family" : "Helvetica", 'vertical-align': 'center', "margin-bottom":"2px"}
extraction_box_style={'padding': '20px','background-color': '#f9f9f9','border-radius': '10px','box-shadow': '0 4px 6px rgba(0, 0, 0, 0.1)','flex': '1','min-width': '250px', 'flex-basis': '300px', 'height':'100px'}
extraction_box_cl_style={'padding': '20px','max-height':'240px','overflow-y':'auto','background-color': '#f9f9f9','border-radius': '10px','box-shadow': '0 4px 6px rgba(0, 0, 0, 0.1)','flex': '1','min-width': '250px', 'flex-basis': '300px', 'height':'290px'}
cool_button_style={'font-size': '24px','padding': '30px 60px','background': 'linear-gradient(135deg, #4CAF50, #81C784)','color': 'white','border': 'none','border-radius': '10px','box-shadow': '0px 4px 6px rgba(0, 0, 0, 0.1)','cursor': 'pointer','transition': 'transform 0.2s, box-shadow 0.2s','flex':'1','min-width': '300px', 'flex-basis': '350px','height': '290px'}
stop_button_style = {'font-size': '18px','padding': '30px 60px','background': 'linear-gradient(135deg, #FF0000, #FF6347)','color': 'white','border': 'none','border-radius': '10px','box-shadow': '0px 4px 6px rgba(0, 0, 0, 0.1)','cursor': 'pointer','transition': 'transform 0.2s, box-shadow 0.2s','flex': '1','min-width': '300px','flex-basis': '350px','height': '200px'}
extraction_cl_style={'margin-top':'50px','max-height':'750px','overflow-y':'auto','background-color': '#f9f9f9','border-radius': '10px','box-shadow': '0 4px 6px rgba(0, 0, 0, 0.1)','flex': '1','min-width': '250px', 'flex-basis': '300px', 'display':'flex', 'align-items':'flex-start'}
# Color definitions
color_points = 'cornflowerblue'

# Directory of the script
SCRIPT_DIR = os.path.dirname(__file__)


# Dash app initialization
# app = dash.Dash(__name__)

import dash_bootstrap_components as dbc

app = dash.Dash(
    __name__,
    external_stylesheets=[dbc.themes.FLATLY],
    meta_tags=[{"name": "viewport", "content": "width=device-width, initial-scale=1"}],
    assets_folder='BFmonitor/assets'
)
app.title = "BinderFlow Monitor"
server = app.server

# Define layout
def serve_layout():
    return html.Div(
        className="app-wrapper",
        children=[
            # Sidebar
            html.Div(
                className="sidebar",
                children=[
                    html.H1("BinderFlow Monitor", className="app-title"),
                    # Removed: html.Div(id='row-count'),
                    html.Hr(),
                    html.H4("Project Directory"),
                    html.Pre(path_to_tree(working_dir), style={
                        "fontFamily": "monospace",
                        "whiteSpace": "pre",
                        "backgroundColor": "#f4f4f4",
                        "padding": "10px",
                        "borderRadius": "6px"
                    }),
                    html.Hr(),
                    html.H4("Create Campaign"),
                    dbc.Row([
                        dbc.Col(
                            dcc.Input(id='campaign-name', type='text', placeholder='Enter campaign name'),
                            className='dropdown', style ={'width':'100%','margin-right':'5px', 'textAlign':'center'}
                        ),
                        dbc.Col(
                            dbc.Button('Create Campaign', id='create-campaign-button', n_clicks=0, className="button"),
                            style={'textAlign':'center', 'margin-top':'10px'}
                        ),
                    ]),
                    html.Hr(),
                    html.H4("Working Campaign Directory:"),
                    dcc.Dropdown(
                        options=[
                            {
                                'label': (
                                    f"{os.path.basename(d)}"
                                ),
                                'value': d
                            }
                            for d in directories_list
                        ],
                        value=directories_list[0],
                        placeholder='Select a folder',
                        id='directory-dropdown',
                        className='dropdown'
                    ),
                    html.Br(),
                    html.Button('STOP CAMPAIGN', id='stop-campaign', n_clicks=0, className="button button-danger"),
                    html.Br(),
                    html.Button('Filters & Axes', id='open-filters', n_clicks=0, className='button'),
                    dbc.Collapse(
                        html.Div([
                            html.H5('PAE Interaction'),
                            dcc.RangeSlider(id='pae_interaction_thres', min=0, max=30, value=[0,10], tooltip={"placement":"bottom","always_visible":True}),
                            html.H5('CUTRE'),
                            dcc.RangeSlider(id='CUTRE_thres', min=0, max=70, value=[0,10], tooltip={"placement":"bottom","always_visible":True}),
                            html.H5('pLDDT binder'),
                            dcc.RangeSlider(id='plddt_binder_thres', min=0, max=100, value=[80,100], tooltip={"placement":"bottom","always_visible":True}),
                            html.H5('dSASA'),
                            dcc.RangeSlider(id='dsasa_thres', min=0, max=10000, value=[1000,10000], tooltip={"placement":"bottom","always_visible":True}),
                            html.H5('Shape Complementarity'),
                            dcc.RangeSlider(id='shape_complementarity_thres', min=0, max=1, value=[0.5,1], tooltip={"placement":"bottom","always_visible":True}),
                            html.H5('Interface HBond'),
                            dcc.RangeSlider(id='interface_hbond_thres', min=0, max=15, value=[3,15], tooltip={"placement":"bottom","always_visible":True}),
                            html.H5('Interface Unsatisfied HBond'),
                            dcc.RangeSlider(id='interface_unsat_hbond_thres', min=0, max=15, value=[0,4], tooltip={"placement":"bottom","always_visible":True}),
                            html.H5('Binder Surface Hydrophobicity'),
                            dcc.RangeSlider(id='binder_surf_hyd_thres', min=0, max=1, value=[0,0.35], tooltip={"placement":"bottom","always_visible":True}),
                            html.H5('ipSAE'),
                            dcc.RangeSlider(id='ipSAE', min=0, max=1, value=[0.6,1], tooltip={"placement":"bottom","always_visible":True}),
                            
                            html.Hr(),
                            html.Div([
                                html.H4('Axis Values')]),
                            html.Div([
                                html.Span('X', className='axis-label'),
                                dcc.Dropdown([
                                                'plddt_binder',
                                                'pae_interaction',
                                                'CUTRE',
                                                'dG',
                                                'dSASA',
                                                'Shape_complementarity',
                                                'Packstat',
                                                'dG_SASA_ratio',
                                                'length',
                                                'SAP',
                                                'binder_int_hyd',
                                                'binder_surf_hyd',
                                                'interface_hbonds',
                                                'interface_unsat_hbonds',
                                                'ipSAE',
                                                'RMSD'
                                            ], 'pae_interaction', id='xaxis_value', className='dropdown', style = {'width': '100%'})
                            ], className='axis-control'),
                            html.Div([
                                html.Span('Y', className='axis-label'),
                                dcc.Dropdown( [
                                                'plddt_binder',
                                                'pae_interaction',
                                                'CUTRE',
                                                'dG',
                                                'dSASA',
                                                'Shape_complementarity',
                                                'Packstat',
                                                'dG_SASA_ratio',
                                                'length',
                                                'SAP',
                                                'binder_int_hyd',
                                                'binder_surf_hyd',
                                                'interface_hbonds',
                                                'interface_unsat_hbonds',
                                                'ipSAE',
                                                'RMSD'
                                            ], 'plddt_binder', id='yaxis_value', className='dropdown', style= {'width': '100%'})
                            ], className='axis-control'),
                        ], style={'padding':'10px'}),
                        id='filters-collapse',
                        is_open=False
                    ),
                    html.Br(),
                ]
            ),
            # Main content
            html.Div(
                className="main-content",
                children=[
                    html.Div([
                        dcc.Tabs([
                            dcc.Tab(label='Input', children=[
                                # molstar, upload for pdb and json, check_config, execute bf and pipeline traking
                                dbc.Row([
                                    dbc.Col([
                                        dbc.Card([
                                            dbc.CardHeader("PDB preview", className='card-header-primary'),
                                            dbc.CardBody([
                                                dcc.Store(id='selected-residues', data='No residues'),
                                                dbc.Row([
                                                    dbc.Col([
                                                        dbc.Col([
                                                            dcc.Store(id='selected-process', data='Protein_Design'),
                                                            dcc.Store(id='pd_dict', data='No_data'),
                                                            dbc.DropdownMenu(label='Select Model',
                                                                        children=[
                                                                            dbc.DropdownMenuItem("Protein Design", id={'type':'input-model', 'name':'Protein_Design'}, n_clicks=0),
                                                                            dbc.DropdownMenuItem("Partial Diffusion", id={'type':'input-model', 'name':'Partial_Diffusion'}, n_clicks=0),
                                                                            dbc.DropdownMenuItem("Sequence Diversity", id={'type':'input-model', 'name':'Sequence_Diversity'}, n_clicks=0)
                                                                        ],
                                                                        id='process-dropdown')
                                                            ]),
                                                        dbc.Row([
                                                            dbc.Col([
                                                        dcc.Upload(
                                                        id='upload-input-pdb',
                                                        children=html.Div([
                                                            'Input ',
                                                            html.A('PDB Files')
                                                        ]),
                                                        style={
                                                            'width': '100%',
                                                            'height': '60px',
                                                            'lineHeight': '60px',
                                                            'borderWidth': '1px',
                                                            'borderStyle': 'dashed',
                                                            'borderRadius': '5px',
                                                            'textAlign': 'center',
                                                            'margin': '10px'
                                                        },
                                                        # Allow multiple files to be uploaded
                                                        multiple=True
                                                        )]),
                                                        dcc.Store(id='selected_template', data='None'),
                                                        dbc.Col([
                                                        dcc.Upload(
                                                        id='upload-template-pdb',
                                                        children=html.Div([
                                                            'Template ',
                                                            html.A('PDB Files')
                                                        ]),
                                                        style={
                                                            'width': '100%',
                                                            'height': '60px',
                                                            'lineHeight': '60px',
                                                            'borderWidth': '1px',
                                                            'borderStyle': 'dashed',
                                                            'borderRadius': '5px',
                                                            'textAlign': 'center',
                                                            'margin': '10px'
                                                        },
                                                        # Allow multiple files to be uploaded
                                                        multiple=True
                                                    )])])
                                                ])
                                                ]),
                                                dashmolstar.MolstarViewer(
                                                                        id="input_pdb_molecule",
                                                                        style={'height': '600px', 'width': '100%'},
                                                                        data=[]
                                                                        ),
                                                dbc.Row([
                                                    dbc.Col([
                                                        dcc.Store(id='selected-surface', data='cartoon'),
                                                        dbc.DropdownMenu(label='Select Model',
                                                                    children=[
                                                                        dbc.DropdownMenuItem("gaussian-surface", id={'type':'surf-model', 'name':'gaussian-surface'}, n_clicks=0),
                                                                        dbc.DropdownMenuItem("cartoon", id={'type':'surf-model', 'name':'cartoon'}, n_clicks=0)
                                                                    ],
                                                                    id='surface-dropdown')
                                                    ]),
                                                ]),
                                                dbc.Row([
                                                    dbc.Col([
                                                    html.Div(
                                                        dbc.Button('Get Selected elements', id='get-input-elements', n_clicks=0, className='button'),
                                                        style={'textAlign':'center', 'margin-top':'10px'}
                                                        )]),
                                                    dbc.Col([
                                                    html.Div(
                                                        dbc.Button('Crop Model with selected elements', id='corp_model', n_clicks=0, className='button'),
                                                        style={'textAlign':'center', 'margin-top':'10px'}
                                                        )])
                                                        ]),
                                                html.Div(id='residue-viewer', children="Waiting for residue selection...", className='outfile_viewer', style={'maxHeight':'185px'}),
                                                html.Hr(),
                                                html.H6('Binder size'),
                                                dcc.RangeSlider(id='binder_size', min=5, max=500, value=[50,150], tooltip={"placement":"bottom","always_visible":True})
                                            ])
                                        ], className='mb-4'),
                                        html.Hr(),
                                    ], width = 8),
                                    dbc.Col([
                                        dbc.Row([
                                            dbc.Card([
                                                dbc.CardHeader("JSON Preview", className='card-header-primary', style={'width':'100%'}),
                                                dbc.CardBody([
                                                    html.Div(id='json-viewer', children="Waiting for JSON...", className='outfile_viewer', style={'maxHeight':'625px'})
                                                ])
                                            ], className='mb-4', style={'height':700, 'width':'100%'})]),
                                        dbc.Row([
                                            dcc.Upload(
                                                    id='upload-input-json',
                                                    children=html.Div([
                                                        'Drag and Drop or ',
                                                        html.A('Select input JSON Files')
                                                    ]),
                                                    style={
                                                        'width': '100%',
                                                        'height': '60px',
                                                        'lineHeight': '60px',
                                                        'borderWidth': '1px',
                                                        'borderStyle': 'dashed',
                                                        'borderRadius': '5px',
                                                        'textAlign': 'center',
                                                        'margin': '10px'
                                                    },
                                                    # Allow multiple files to be uploaded
                                                    multiple=True
                                                ),
                                        ]),
                                        dbc.Row([
                                            dbc.Row([
                                                dbc.Col([
                                                    html.H5('Max Threads'),
                                                    dcc.Input(id='max_threads',type='number',placeholder='Maximun number of nodes used in parallel',className="input-box", style={'width':'100%'}, min=1, max=4),
                                                ]),
                                                dbc.Col([
                                                    html.H5('RFD Designs'),
                                                    dcc.Input(id='rfd_designs',type='number',placeholder='Number of backbones to design',className="input-box", style={'width':'100%'}, min=1, max=100),
                                                ]),
                                            ]),
                                            dbc.Row([
                                                dbc.Col([
                                                html.H5('pMPNN nseqs'),
                                                dcc.Input(id='pmpnn_seqs',type='number',placeholder='Number of sequences to design',className="input-box", style={'width':'100%'}, min=1),
                                                ]),
                                                dbc.Col([
                                                html.H5('pMPNN relax cycles'),
                                                dcc.Input(id='pmpnn_cycles',type='number',placeholder='Number of cycles (i.e. 0)',className="input-box", style={'width':'100%'}, min=0, max=1),
                                                ]),
                                            ]),
                                            dbc.Row([
                                                dbc.Col([
                                                html.H5('Noise steps'),
                                                dcc.Input(id='noise_steps',type='number',placeholder='Number of steps of noise added',className="input-box", style={'width':'100%'}, min=10),
                                                ]),
                                                dbc.Col([
                                                html.H5('Noise scale'),
                                                dcc.Input(id='noise_scale',type='text',placeholder='Amount of noise added at each step',className="input-box", style={'width':'100%'},min=0.1, max=2, step=0.1),
                                                ]),
                                            ]),
                                            dbc.Row([
                                                dbc.Col([
                                                html.H5('Checkpoint'),
                                                dcc.Dropdown(id='ckp',options=[],className='dropdown', style={'width':'100%'}),
                                                ]),
                                                dbc.Col([
                                                html.H5('Node'),
                                                dcc.Input(id='node',type='text',placeholder='',className="input-box", style={'width':'100%'}),
                                                ]),
                                            ]),
                                            dbc.Row([
                                                dbc.Col([
                                                    html.H5('Residues'),
                                                    dcc.Input(id='residue',type='text',placeholder='Residues to fixed in the structure',className="input-box", style={'width':'100%'}),
                                                ]),
                                                dbc.Col([
                                                    html.H5('Hits number'),
                                                    dcc.Input(id='hits',type='number',placeholder='Nº of hits desired',className="input-box", style={'width':'100%'}, min=1, max=1000),
                                                ]),
                                            ])
                                        ]),
                                        dbc.Row([
                                            html.Div(
                                                dbc.Button('Generate JSON', id='generate-json-button', n_clicks=0, className='button'),
                                                style={'textAlign':'center', 'margin-top':'10px'}
                                            ) 
                                        ]),

                                    ]),
                                ]),

                                html.Hr(),
                                dbc.Row([
                                    dbc.Col([
                                        dbc.Card([
                                        dbc.CardHeader("Binderflow Initialize", className='card-header-af3'),
                                        dbc.CardBody([
                                            html.Div(id='binderflow-initial-viewer', children="Waiting to execute Binderflow...", className='outfile_viewer', style={'maxHeight':'185px'}),
                                        html.Div(
                                            dbc.Button('Execute Binderflow', id='execute-binderflow-button', n_clicks=0, className='button_bindeflow'),
                                            style={'textAlign':'center', 'margin-top':'10px'}
                                        )])
                                    ], className='mb-4',style={'height':'100%', 'width':'100%'}),
                                    ])
                                ])
                            ]),
                            dcc.Tab(label='Live Watcher', children=[
                                # Metrics cards row
                                html.Div(
                                    children=[
                                        dbc.Card(
                                            [
                                                html.Div(id='row-count', className='metric-value'),
                                                html.Div("Finished Designs", className='metric-label'),
                                            ],
                                            className="metric-card"
                                        ),
                                        dbc.Card(
                                            [
                                                html.Div(id='hit-count', className='metric-value'),
                                                html.Div("Hits", className='metric-label'),
                                            ],
                                            className="metric-card"
                                        ),
                                        dbc.Card(
                                            [
                                                html.Div(id='hit-efficiency', className='metric-value'),
                                                html.Div("Hit Efficiency", className='metric-label'),
                                            ],
                                            className="metric-card"
                                        ),
                                    ],
                                    className="metrics-row"
                                ),
                                html.Div([
                                    html.Div([
                                        html.Div([
                                            # Removed filter/axis toggles and panels; controls now in Offcanvas
                                            dbc.Card(
                                                [
                                                    dbc.CardHeader("Scatter Plot"),
                                                    dbc.CardBody(
                                                        dcc.Graph(
                                                            id='scatter-plot',
                                                            className='graph-container',
                                                            config={'toImageButtonOptions': {'format': 'svg', 'filename': 'scatter_plot', 'height': 600, 'width': 800, 'scale': 1}, 'displaylogo': False}
                                                        )
                                                    ),
                                                ],
                                                className="graph-card"
                                            ),
                                            dbc.Card(
                                                [
                                                    dbc.CardHeader("Radar Plot"),
                                                    dbc.CardBody(
                                                        dcc.Graph(
                                                            id='radar-plot',
                                                            className='graph-container',
                                                            config={'toImageButtonOptions': {'format': 'svg', 'filename': 'radar_plot', 'height': 600, 'width': 800, 'scale': 1}, 'displaylogo': False}
                                                        )
                                                    ),
                                                ],
                                                className="graph-card"
                                            ),
                                        ], className='plot-and-controls'),
                                    ], className='plot-and-controls'),
                                ]),
                                html.Div(
                                    [
                                        html.Div(
                                            [
                                                dcc.Store(id='input_pdb_path', data='No_path'),
                                                dcc.Upload(
                                                    id='upload-pd',
                                                    children=html.Div([
                                                        'Drag and Drop or ',
                                                        html.A('Select input PDB Files for comparasion')
                                                    ]),
                                                    style={
                                                        'width': '100%',
                                                        'height': '60px',
                                                        'lineHeight': '60px',
                                                        'borderWidth': '1px',
                                                        'borderStyle': 'dashed',
                                                        'borderRadius': '5px',
                                                        'textAlign': 'center',
                                                        'margin': '10px'
                                                    },
                                                    # Allow multiple files to be uploaded
                                                    multiple=True
                                                )
                                            ],
                                            style={'width':'100%', 'padding':'10px'}
                                        )
                                    ],
                                    style={'display': 'flex', 'flex-direction': 'row', 'align-items': 'center'}
                                ),
                            ]),
                            dcc.Tab(label='Pipeline Tracking', children=[
                                dbc.Row([
                                    dbc.Col(
                                        dbc.Card([
                                            dbc.CardHeader(id='job-status-counts', className='card-header-primary'),
                                            dbc.ListGroup(id='job-status-list', flush=True, className='mb-4', style={'width':'100%', 'overflowY': 'auto', 'maxHeight': '750px'}),
                                            dcc.Interval(id='interval-component', interval=60000, n_intervals=0),
                                            dcc.Interval(id="interval-binderflow", interval=3000, n_intervals=0)],className='outfile_cards'
                                            ),
                                            width=4
                                    ),
                                    dbc.Col([
                                        dbc.Row([
                                            dbc.Col(
                                                dbc.Card([
                                                    dbc.CardHeader("Outfile Viewer", className='card-header-primary'),
                                                    dbc.CardBody([
                                                        html.Div(id='log-viewer', className='outfile_viewer', style={'maxHeight':'325px'}),
                                                        html.Div(id='log-viewer-content', style={'display': 'none'})
                                                ]),
                                                ], className='outfile_card',style={'height': 400})
                                            )                                      
                                        ]),
                                        dbc.Row([
                                            dbc.Col(
                                                dbc.Card([
                                                    dbc.CardHeader("Errorfile Viewer", className='card-header-primary'),
                                                    dbc.CardBody([
                                                        html.Div(id='error-viewer',className='outfile_viewer', style={'maxHeight':'325px'}),
                                                        html.Div(id='error-viewer-content', style={'display': 'none'})
                                                    ]),
                                                ], className='outfile_card',style={'height': 400})
                                            )
                                        ]),
                                    ], width=8)
                                ]),
                                dbc.Row([
                                    dbc.Col(
                                        dbc.Card([
                                            dbc.CardHeader("Binderflow Viewer", className='card-header-af3'),
                                            dbc.CardBody([
                                                html.Div(id='binderflow-viewer', children="Waiting for project.log...",className='outfile_viewer', style={'maxHeight':'325px'}),
                                        ]),
                                        ], className='outfile_card',style={'height': 400})
                                    ) 
                                ])
                            ]), 
                            dcc.Tab(label='Extraction', children=[
                                # Hit Viewer & Selection
                                dbc.Row([
                                    dbc.Col(
                                        dbc.Card([
                                            dbc.CardHeader("Hit PDB preview", className='card-header-primary'),
                                            dbc.CardBody([
                                                dcc.Dropdown(
                                                    options=[],
                                                    value=None,
                                                    placeholder='Select a hit',
                                                    id='extractions_molecule_dropdown',
                                                    className='dropdown'
                                                ),
                                                dashmolstar.MolstarViewer(
                                                                        id="extractions_molecule",
                                                                        style={'height': '600px', 'width': '100%'},
                                                                        )
                                            ])
                                        ], className='mb-4'),
                                        width=7
                                    ),
                                    dbc.Col(
                                        dbc.Card([
                                            dbc.CardHeader("Selection of hits for extraction", className='card-header-primary'),
                                            dbc.CardBody([
                                                dcc.Checklist(
                                                    options=[],
                                                    value=[],
                                                    id='extraction-selection',
                                                    style={'display':'flex','flexDirection':'column','gap':'5px', 'width':'100%','maxHeight':'600px', 'overflowY':'auto'},
                                                )
                                            ])
                                        ], className='mb-4'),
                                        width=5
                                    )
                                ]),
                                # Extraction Options Card
                                dbc.Card([
                                    dbc.CardHeader("Extraction Model", className='card-header-primary'),

                                    dbc.CardBody([
                                        dbc.Row([
                                            dbc.Col([
                                                dcc.Store(id='selected-model', data='PDB'),
                                                dbc.DropdownMenu(label='Select Model',
                                                                    children=[
                                                                        dbc.DropdownMenuItem("PDB", id={'type':'dna-model', 'name':'PDB'}, n_clicks=0),
                                                                        dbc.DropdownMenuItem("CodonTransformer", id={'type':'dna-model', 'name':'CT'}, n_clicks=0),
                                                                        dbc.DropdownMenuItem("AlphaFold3", id={'type':'dna-model', 'name':'AF3'}, n_clicks=0),
                                                                        dbc.DropdownMenuItem("JohnBercow", id={'type':'dna-model', 'name':'JB'}, n_clicks=0),
                                                                    ],
                                                                    id='dna-model-dropdown'),
                                                html.Div( id='model-specific-options'), # Division to put model-specific options
                                            ], width=6),
                                            dbc.Col([
                                                dbc.Card([
                                                    dbc.CardHeader("Extraction outfile"),
                                                    dbc.CardBody([
                                                        html.Div(id='extraction-viewer', children='Initial content...',style={"maxHeight": "1150px", "overflow": "scroll"}),
                                                    ])
                                                ], className='outfile_viewer', style={'height':'100%'})
                                            ], width=6)
                                        ]),
                                        html.Div(
                                            dbc.Button('Execute operation', id='execute-hits', n_clicks=0, className='button'),
                                            style={'textAlign':'center', 'margin-top':'10px'}
                                        ),
                                        dcc.Interval(
                                            id='interval-component-extraction',
                                            interval=5000,  # Update every 5000 milliseconds (3 second)
                                            n_intervals=0
                                        ),
                                    ])
                                ], className='mb-4'),
    
                            ]),
                        ])
                    ])
                ]
            ),
            # Offcanvas removed
        ]
    )

app.layout = serve_layout

#########################
## INPUT TAB CALLBACKS ##
#########################

# callback to create campaigns directories
@callback(
    Output('directory-dropdown', 'options'),
    Input("create-campaign-button", "n_clicks"),
    State("campaign-name", "value"),
    prevent_initial_call=True
)
def create_campaign(n_clicks, campaign_name):
    ctx = dash.callback_context
    if ctx.triggered_id == "create-campaign-button":
        create_campaigns(working_dir, campaign_name)
        directories_list = get_working_directories(working_dir)
        options = [{'label': os.path.basename(dir), 'value': dir} for dir in directories_list]
        return options

#callback for design
@callback(
    Output('selected-process', 'data'),
    Output('process-dropdown', 'label'),
    Input({'type':'input-model', 'name': dash.ALL}, 'n_clicks'),
)
def update_process(n):
    ctx = dash.callback_context
    if not ctx.triggered:
        selected_model='Protein_Design'
        return ('Protein_Design', 'Protein Design')
    else:
        triggered_id = ctx.triggered[0]['prop_id'].split('.')[0]
        selected_model = eval(triggered_id)['name']

    if selected_model == "Protein_Design":
        return ('Protein_Design', 'Protein Design')
    elif selected_model == "Partial_Diffusion":
        return ('Partial_Diffusion', 'Partial Diffusion')
    elif selected_model == "Sequence_Diversity":
        return ('Sequence_Diversity', 'Sequence Diversity')

#callback to execute binderflow
@callback(
    Output("binderflow-initial-viewer","children"),
    Input("execute-binderflow-button","n_clicks"),
    State('directory-dropdown', 'value'),
    State('selected-process', "data"),
    State("pd_dict", "data"),
    prevent_initial_call=True    
)
def perform_binderflow(click, dir,process, pd_dict_str):
    import time
    ctx = dash.callback_context
    if ctx.triggered_id == "execute-binderflow-button":
        print('Running Binderflow in directory '+str(dir)+' with process '+str(process))
        message = run_binderflow(SCRIPT_DIR, dir, pd_dict_str, process)
        return message
        

#callback for project reader
@callback(
        Output("binderflow-viewer","children"),
        Input("interval-binderflow","n_intervals"),
        State('directory-dropdown', 'value'),
        State('selected-process', "data"),
)
def update_binderflow_viewer(n,dir, process):
    try:
        with open(str(dir)+"/project_binderflow_Design.log","r") as file:
            return file.read()
    except:
        return str("Waiting for log file...")
    
#callback for surface
@callback(
    Output('selected-surface', 'data'),
    Output('surface-dropdown', 'label'),
    Input({'type':'surf-model', 'name': dash.ALL}, 'n_clicks'),
)
def update_surface(n):
    ctx = dash.callback_context
    if not ctx.triggered:
        selected_model='cartoon'
        return ('cartoon', 'cartoon')
    else:
        triggered_id = ctx.triggered[0]['prop_id'].split('.')[0]
        selected_model = eval(triggered_id)['name']

    if selected_model == "gaussian-surface":
        return ('gaussian-surface', 'gaussian-surface')
    elif selected_model == "cartoon":
        return ('cartoon', 'cartoon')
    
#callback for input pdb
@callback(
        Output("input_pdb_molecule","data"),
        Output('upload-input-pdb', 'contents'),
        Output("pd_dict", "data"),
        Input("upload-input-pdb","contents"),
        Input({'type':'surf-model', 'name': dash.ALL}, 'n_clicks'),
        Input("corp_model","n_clicks"),
        State("input_pdb_molecule","data"),
        State('upload-input-pdb', 'filename'),
        State('directory-dropdown', 'value'),
        State('selected-surface', "data"),
        State('selected-process', "data"),
        State("input_pdb_molecule","selection"),
        prevent_initial_call=True
)
 
def input_pdb(pdb_file, surface,clk,stored_data, filename,dir,surf, process, residues):
    ctx = dash.callback_context
    triggered = ctx.triggered_id
    pd_dict_string="No_data"
    if triggered == 'upload-input-pdb':
        data = load_input_pdb(filename, process, pdb_file, dir, surf)
        
        if process == "Partial_Diffusion" and str(filename[0]).endswith(".pdb"):
            decoded = decoded_pdbFILE(pdb_file)
            pd_dict={}
            control=0
            for line in decoded.split("\n"):
                elements=line.split(" ")
                if elements[0] == "CUTRE":
                    control=1
                    pd_dict[elements[0]] = elements[1]
                elif control == 1 and line != "" and line != "\n":
                    pd_dict[elements[0]] = elements[1]
            pd_dict_string=json.dumps(pd_dict)
            print(pd_dict_string)
        return data, None, pd_dict_string
    
    elif triggered == "corp_model" and process == "Protein_Design":
        target = Target(residues)
        sout="["
        #chains=len(target)
        for ch in target.chains:
            for res in ch.residues:
                sout=sout+str(ch.name)+str(res.number)+","
        if sout != "[":
            sout=sout[:-1]+"]"
        else: 
            sout="No residues"
        if sout != "No residues":
            corp_PDB(sout, dir)
            chainB = molstar_helper.get_targets(chain="B")
            chainB_representation = Representation(type="cartoon", color="hydrophobicity")
            component_B = molstar_helper.create_component("Chain B", chainB, chainB_representation)
            components=[component_B]

            # Load the molecule data
            pdb_file=f"{dir}/input/Protein_Design_input.pdb"
            data = molstar_helper.parse_molecule(pdb_file, fmt='pdb',component=components)

            return [data], None, "No_data"
        else:
            return dash.no_update, dash.no_update, dash.no_update
    # user changed surface model
    if isinstance(triggered, dict) and triggered.get('type') == 'surf-model':
        target_surface = json.loads(ctx.triggered[0]['prop_id'].split(".")[0])['name'] # This line overcomes problems with the number of clicks
        print(f"Surface changed to: {target_surface} \n")
        # Modify stored_data representation key to the new surface
        stored_data[0]['component'][0]['representation'] = [{'type': target_surface, 'color': 'hydrophobicity'}]
        data = stored_data
        return data, dash.no_update, dash.no_update

    return dash.no_update, dash.no_update, dash.no_update


#callback for input json
@app.callback(Output('json-viewer', 'children'),
            Output('upload-input-json', 'contents'),
            Input('generate-json-button',"n_clicks"),
            Input('upload-input-json', 'contents'),
            State('upload-input-json', 'filename'),
            State('selected-process', "data"),
            State('selected-residues', 'data'),
            State("max_threads","value"),
            State('rfd_designs',"value"),
            State('pmpnn_seqs',"value"),
            State('pmpnn_cycles',"value"),
            State('noise_steps',"value"),
            State('noise_scale',"value"),
            State("ckp", "value"),
            State("node", "value"),
            State('residue', 'value'),
            State('hits', "value"),
            State('directory-dropdown', 'value'),
            State("binder_size","value"),
            State("selected_template","data"),
            prevent_initial_call=True
)
def update_input_json(n, contents, list_of_names,process, residues_molstar, max_threads, rfd_designs,pmpnn_seqs,pmpnn_cycles, noise_steps,noise_scale,ckp,node,residue,hits,dir, size, template):
    ctx = dash.callback_context
    SCRIPT_DIR= os.path.dirname(__file__)
    if contents is not None:
        try:
            json_string=json_uploader(contents,dir)
            return (json_string, None)
        except:
            return (str(str(list_of_names)+" is not a valid JSON file"), None)
    elif ctx.triggered_id == "generate-json-button":
        json_file=generate_json_from_flags(process, size, dir, template, max_threads, residues_molstar, rfd_designs, pmpnn_seqs, pmpnn_cycles, noise_steps, noise_scale, ckp, hits, SCRIPT_DIR)
        return (json_file, None)

#callback for template pdb
@callback(
        Output("selected_template","data"),
        Output('upload-template-pdb', 'contents'),
        Input("upload-template-pdb","contents"),
        State('upload-template-pdb', 'filename'),
        State('directory-dropdown', 'value'),
        State('selected-process', "data"),
        prevent_initial_call=True
)
def get_template(template, filename, dir, process):
    if template is not None and process != "Partial_Diffusion":
        try: 
            pdb_path = load_template_pdb(process,template,dir)
            process_pdb_B_1000(pdb_path, dir, process, filename, template=True)
            return (f'{dir}/input/{process}_template.pdb', None)
        except:
            return (str(str(filename)+" is not a valid PDB file"), None)

#callback to save selected residues
@callback(
        Output("selected-residues","data"),
        Output("residue-viewer","children"),
        Input("get-input-elements","n_clicks"),
        State("input_pdb_molecule","selection"),
        State('directory-dropdown', 'value'),
        State('selected-process', "data"),
        prevent_initial_call=True
)
def get_residues(clk, residues,dir, process):
    ctx = dash.callback_context
    if ctx.triggered_id == "get-input-elements" and process == "Protein_Design":
        target = Target(residues)
        sout="["
        #chains=len(target)
        with open(str(dir)+"/pdb_selection.json","w") as f:
            for ch in target.chains:
                f.write("Chain "+str(ch.name)+": ")
                for res in ch.residues:
                    f.write(str(res.name)+str(res.number)+", ")
                    sout=sout+str(ch.name)+str(res.number)+","
                f.write("\n")
        if sout != "[":
            sout=sout[:-1]+"]"
        else: 
            sout="No residues"
        return (sout, sout)
    else:
        return ("No residues", "No residue selection for Partial Diffusion")
        
# Function to get the checkpoint paths

@callback(
    Output('ckp', 'options'),
    Input('interval-component', 'n_intervals'),
)
def get_checkpoint_paths(intervals):
    checkpoint_files, _ = get_ckp_paths()
    return checkpoint_files



####################################
#### LIVE WATCHER TAB CALLBACKS ####
####################################

# Callback to update graphs and table
@callback(
    Output('scatter-plot', 'figure'),
    Output('row-count', 'children'),
    Output('hit-count', 'children'),
    Output('hit-efficiency', 'children'),
    Output('extractions_molecule_dropdown', 'options'),
    Output('extraction-selection', 'options'),
    Output('extraction-selection', 'value'),
    Output('input_pdb_path', 'data'),
    Output('upload-pd', 'contents'),
    Input("upload-pd","contents"),
    Input('upload-pd', 'filename'),
    Input('directory-dropdown', 'value'),
    Input('interval-component', 'n_intervals'),
    Input('xaxis_value', 'value'),
    Input('yaxis_value','value'),
    Input('pae_interaction_thres', 'value'),
    Input('CUTRE_thres', 'value'),
    Input('plddt_binder_thres', 'value'),
    Input('dsasa_thres', 'value'),
    Input('shape_complementarity_thres', 'value'),
    Input('interface_hbond_thres', 'value'),
    Input('interface_unsat_hbond_thres', 'value'),
    Input('binder_surf_hyd_thres', 'value'),
    Input('ipSAE','value'),
    
)

def update_graph(contents, filename,working_dir, n, xaxis_value, yaxis_value, pae_interaction_thres, CUTRE_thres, plddt_binder_thres, dsasa_thres, shape_complementarity_thres, interface_hbond_thres, interface_unsat_hbond_thres, binder_surf_hyd_thres,ipsae):
    ctx = dash.callback_context
    directory = f'{working_dir}/output/'
    input_pdb_path=None
    # adding input pdb path as popup
    if contents is not None and str(filename[0]).endswith(".pdb"):
            #name=filename[0]
            pattern = os.path.join(working_dir, '**', filename[0])
            pdb_path = glob.glob(pattern, recursive=True)[0]
            name = pdb_path.replace(working_dir,"")[1:]
            print(name)
            if os.path.isfile(f"{working_dir}/{name}") == True:
                input_pdb_path=name
            else:
                print("No path found")
    
    merged_df = pd.DataFrame()
    merged_df = merge_csv_files(working_dir, input_pdb_path)
    if merged_df.empty:
        scatter_plot = go.Figure()
        scatter_plot.add_trace(go.Scatter(x=[], y=[], mode='markers'))
        return scatter_plot, "No designs finished yet", "0", "0%", [], [], []
    filtered_df = filtering_df(merged_df, pae_interaction_thres, CUTRE_thres, plddt_binder_thres, dsasa_thres, shape_complementarity_thres, interface_hbond_thres, interface_unsat_hbond_thres, binder_surf_hyd_thres,ipsae)


    # Make scatter plot
    scatter_plot, row_count_text = update_scatter_plot(working_dir, merged_df, filtered_df, xaxis_value, yaxis_value, input_pdb_path)



    # Dropdown HITS
    dropdown_options = get_hit_names(filtered_df, xaxis_value)
    
    # Metrics for cards
    finished_models = len(merged_df)
    hit_count = len(dropdown_options)
    if dropdown_options == ["No hits found using current filters"]:
        hit_count = 0
        hit_efficiency = "0%"
    elif finished_models > 0:
        hit_efficiency = f"{(hit_count/finished_models*100):.1f}%"
    else:
        hit_efficiency = "N/A"

    return scatter_plot, row_count_text, hit_count, hit_efficiency, dropdown_options, dropdown_options, dropdown_options, input_pdb_path, None
    


#callback to update the radar plot
@callback(
    Output('radar-plot', 'figure'),
    [
    Input('scatter-plot', 'clickData'),
    Input('interval-component', 'n_intervals'),
    Input('directory-dropdown', 'value'),
    Input('input_pdb_path', 'value'),
    Input('pae_interaction_thres', 'value'),
    Input('CUTRE_thres', 'value'),
    Input('plddt_binder_thres', 'value'),
    Input('dsasa_thres', 'value'),
    Input('shape_complementarity_thres', 'value'),
    Input('interface_hbond_thres', 'value'),
    Input('interface_unsat_hbond_thres', 'value'),
    Input('binder_surf_hyd_thres', 'value'),
    ])

#This prints the radar plot
def update_radar_plot(design_to_plot, n, directory, input_pdb_path, pae_interaction_thres,CUTRE_thres, plddt_binder_thres, dsasa_thres, shape_complementarity_thres, interface_hbond_thres, interface_unsat_hbond_thres, binder_surf_hyd_thres):
    #This is meant to store the original designs for the following plotting; is a little bit cutre 
    update_designs_list(designs_list, design_to_plot)
    merged_df = pd.DataFrame()
    merged_df = merge_csv_files(directory, input_pdb_path)
    #Most of the heavy work is carry in the utils file, go there for further info
    radar_figure=radar_plot(designs_list, merged_df,pae_interaction_thres,CUTRE_thres, plddt_binder_thres, dsasa_thres, shape_complementarity_thres, interface_hbond_thres, interface_unsat_hbond_thres, binder_surf_hyd_thres )
    # Override radar area colors with transparency
    return radar_figure


#########################################
#### PIPELINE TRACKING TAB CALLBACKS ####
#########################################


# callback for job-status-list
@callback(
    Output('job-status-list', 'children'),
    Input('interval-component', 'n_intervals'),
    Input('directory-dropdown', 'value')
)
def update_job_list(n, working_dir):
    color = {
        'Waiting':'#C5CBD3',
        'RFD': '#4b2362',
        'Filtering': '#863071',
        'pMPNN':'#c14168',
        'Scoring':'#e5715e',
        'Finished':'#edb081',
        'Failed':'#000F08' 
    }

    df = track_job_status(f'{working_dir}/output/')
    # Extract run number and sort descending
    df_sorted = df.copy()
    if df_sorted.empty:
        return [dbc.ListGroupItem("This is an old project. No job status information available.")]
    df_sorted['run_number'] = df_sorted['job'].str.extract(r'(\d+)', expand=False).astype(int)
    df_sorted = df_sorted.sort_values(['run_number','gpu'], ascending=[False,True ])
    items = []
    for job in df_sorted['job'].unique():
        items.append(
            dbc.ListGroupItem(
                dbc.Row([
                    dbc.Col(
                        html.Span(job, className='me-2', style={'fontWeight': '600'}),
                    width=4),
                    dbc.Col(
                            [
                                dbc.Row(
                                    dbc.DropdownMenu(
                                        label=f'GPU {gpu}: \t {df_sorted["status"][(df_sorted["job"] == job) & (df_sorted["gpu"] == gpu)].values[0]}',
                                        children=[
                                            dbc.DropdownMenuItem('Backbone',  id={'type':'dropdown-item', 'value':f"{df_sorted['path'][(df_sorted['job'] == job) & (df_sorted['gpu'] == gpu)].values[0]}/rfd"}, n_clicks=0),
                                            dbc.DropdownMenuItem('Filtering', id={'type':'dropdown-item', 'value':f"{df_sorted['path'][(df_sorted['job'] == job) & (df_sorted['gpu'] == gpu)].values[0]}/aligning_filtering"}, n_clicks=0),
                                            dbc.DropdownMenuItem('Sequence',  id={'type':'dropdown-item', 'value':f"{df_sorted['path'][(df_sorted['job'] == job) & (df_sorted['gpu'] == gpu)].values[0]}/pmpnn"}, n_clicks=0),
                                            dbc.DropdownMenuItem('Scoring',   id={'type':'dropdown-item', 'value':f"{df_sorted['path'][(df_sorted['job'] == job) & (df_sorted['gpu'] == gpu)].values[0]}/scoring"}, n_clicks=0),
                                        ], color=color[df_sorted['status'][(df_sorted['job'] == job) & (df_sorted['gpu'] == gpu)].values[0]],
                                    )
                                )
                                for gpu in df_sorted['gpu'][df_sorted['job'] == job].unique()
                            ])
                ])
            )
        )
    return items
#Callback to update the log and error viewer
@callback(
    Output('log-viewer', 'children'),
    Output('error-viewer', 'children'),
    Input({"type":'dropdown-item', "value":dash.ALL}, 'n_clicks'),
    prevent_initial_call=True
)

def update_log_viewer(n_clicks_list):
    ctx = dash.callback_context
    if not ctx.triggered:
        return " ", " "
    triggered_id = ctx.triggered[0]['prop_id'].split('.')[0]
    value= eval(triggered_id)['value']

    outfile_path = f"{value}.out"
    errorfile_path = f"{value}.err"
    try:
        with open(outfile_path, 'r') as outfile:
            outfile_content = outfile.read()
    except FileNotFoundError:
        outfile_content = f'{outfile_path} cannot be accessed'
    try:
        with open(errorfile_path, 'r') as errorfile:
            errorfile_content = errorfile.read()
    except FileNotFoundError:
        errorfile_content = f'{errorfile_path} cannot be accessed'

    return outfile_content, errorfile_content

##################################
#### EXTRACTION TAB CALLBACKS ####
##################################

# Callback for molecule representation
@callback(
    Output("extractions_molecule", "data"),
    [
    Input("extractions_molecule_dropdown", "value"),
    Input("directory-dropdown", "value"),
    Input("input_pdb_path", "value")
    ]
)
def return_molecule(value, directory, input_path):
    if not value or not directory:
        raise PreventUpdate

    data_path, filename = get_design_file_path_and_name(value, directory, input_path)


    file_path = os.path.join(data_path, filename + ".pdb")


    print(f"Loading molecule from: {file_path}")

    chainA = molstar_helper.get_targets(chain="A")
    chainB = molstar_helper.get_targets(chain="B")

    chainA_representation = Representation(type='cartoon', color="uniform")
    chainA_representation.set_color_params({'value': 15577217})
    chainB_representation = Representation(type='gaussian-surface', color="uniform")
    chainB_representation.set_type_params({'alpha': 1})
    chainB_representation.set_color_params({'value': 8815233})

    component_A = molstar_helper.create_component("Chain A", chainA, chainA_representation)
    component_B = molstar_helper.create_component("Chain B", chainB, chainB_representation)

    data = molstar_helper.parse_molecule(file_path, "pdb", component=[component_A, component_B], preset={'kind': 'empty'})

    return data

@callback(
    Output('model-specific-options', 'children'),
    Output('selected-model', 'data'),
    Output('dna-model-dropdown', 'label'),
    Input({'type':'dna-model', 'name': dash.ALL}, 'n_clicks'),
)

# Model specific options
# This is a little bit of a hell, probably we should think on moving it into other script and then exporting it
def update_model_specific_options(n_clicks_list):
    ctx = dash.callback_context
    if not ctx.triggered:
        selected_model='PDB'
        return (html.Div(), 'PDB', 'PDB')
    else:
        triggered_id = ctx.triggered[0]['prop_id'].split('.')[0]
        selected_model = eval(triggered_id)['name']

    if selected_model == "AF3":
        return (
            dbc.Card([
                dbc.CardHeader("AlphaFold3 panel", className='card-header-af3'),

                dbc.CardBody([
                    dbc.Row([
                        dbc.Col([
                            dbc.Card(
                                [
                                    dbc.CardHeader("AF3 Input"),
                                    dbc.CardBody(
                                        [
                                            html.Div(id='af3-viewer', children='Initial content...',style={"maxHeight": "325px", "overflow": "scroll"}) # The card content
                                        ]
                                    )
                                ],
                                style={"height":"100%", "width":"100%"}
                            ),
                        ],width=12)]),
                        html.Br(),
                    dbc.Row([
                        dbc.Col([
                            dbc.Card([
                                dbc.CardHeader("AF3 Scatterplot", className='card-header-af3'),
                                dbc.CardBody([
                                    dbc.Row([
                                        dbc.Col([
                                            dcc.Dropdown(
                                                options=['pae_interaction', 'plddt_binder','plddt_binder_AF3', 'ipSAE_AF3', 'ipSAE'],
                                                value=None,
                                                placeholder='X Axis',
                                                id='af3_x_axis',
                                                className='dropdown'
                                            ),
                                        ]),
                                        dbc.Col([
                                                dcc.Dropdown(
                                                    options=['pae_interaction', 'plddt_binder','plddt_binder_AF3', 'ipSAE_AF3', 'ipSAE'],
                                                    value=None,
                                                    placeholder='Y Axis',
                                                    id='af3_y_axis',
                                                    className='dropdown'
                                                ),
                                        ])
                                    ]),
                                    dbc.Row([
                                        dcc.Graph(
                                            id='AF3-scatterplot',
                                            className='graph-container',
                                        )
                                    ])
                                ])
                            ]),
                    ])]),
                    dcc.Interval(
                        id='interval-component-af3',
                        interval=3000,  # Update every 3000 milliseconds (3 second)
                        n_intervals=0
                    ),
                ])
            ], className='af3'),"AF3", "AlphaFold3")

    elif selected_model == 'CT':
        return (html.Div([
                html.Br(),
                html.Div([
                    html.Span("CodonTransformer Options", className='model-options-title')
                ], style={'textAlign':'center'}),
                html.Br(),
                dbc.Row([
                    dbc.Col(
                        dbc.Card([
                            dbc.CardHeader("Add Initial Methionine"),
                            dbc.CardBody(
                                dbc.RadioItems(
                                    [
                                        {'label': 'True', 'value': True},
                                        {'label': 'False', 'value': False}
                                    ],
                                    value=False,
                                    id={'type':'param', 'name':'add-met'},
                                    inputStyle={'margin-right':'5px'}
                                )
                            )
                        ], className='mb-3'),
                        width=6
                    ),
                    dbc.Col(
                        dbc.Card([
                            dbc.CardHeader("Organism"),
                            dbc.CardBody(
                                dcc.Dropdown(
                                    initial_organisms,
                                    'Escherichia coli general',
                                    id={'type':'param', 'name':'organism'},
                                    className='dropdown'
                                ), style={'width': '100%'}
                            )
                        ], className='mb-3'),
                        width=6
                    )
                ]),
                dbc.Row([
                    dbc.Col(
                        dbc.Card([
                            dbc.CardHeader("3' Overhang"),
                            dbc.CardBody(
                                dcc.Input(
                                    placeholder='Sequence',
                                    value='',
                                    id={'type':'param', 'name':'three_prime_overhang'},
                                    className='input-box',
                                    style={'width': '100%'}
                                )
                            )
                        ], className='mb-3'),
                        width=6
                    ),
                    dbc.Col(
                        dbc.Card([
                            dbc.CardHeader("5' Overhang"),
                            dbc.CardBody(
                                dcc.Input(
                                    placeholder='Sequence',
                                    value='',
                                    id={'type':'param', 'name':'five_prime_overhang'},
                                    className='input-box',
                                    style = {'width': '100%'}
                                )
                            )
                        ], className='mb-3'),
                        width=6
                    )
                ]),
                dbc.Row([
                    dbc.Col(
                        dbc.Card([
                            dbc.CardHeader("Minimum Sequence Length"),
                            dbc.CardBody(
                                dcc.Input(
                                    type='number',
                                    value=300,
                                    id={'type':'param', 'name':'random_sequence'},
                                    className='input-box',
                                    style={'width':'100%'}
                                )
                            )
                        ], className='mb-3'),
                        width=6
                    ),
                    dbc.Col(
                        dbc.Card([
                            dbc.CardHeader("GC % of Random Sequence"),
                            dbc.CardBody(
                                dcc.Input(
                                    type='number',
                                    value=50,
                                    id={'type':'param', 'name':'GC_content'},
                                    className='input-box',
                                    style={'width':'100%'}
                                )
                            )
                        ], className='mb-3'),
                        width=6
                    )
                ]),
                dbc.Row([
                    dbc.Col(
                        dbc.Card([
                            dbc.CardHeader("Check Restriction Sites"),
                            dbc.CardBody(
                                dcc.Checklist(
                                    options=[{'label': enz, 'value': enz} for enz in [
                                        'EcoRI','BamHI','HindIII','NotI','XhoI','PstI','SacI','KpnI',
                                        'SmaI','XbaI','SpeI','NcoI','SalI','ApaI','HaeIII','AluI',
                                        'TaqI','BglII','ClaI','MluI','BsaI'
                                    ]],
                                    value=[],
                                    id={'type':'param', 'name':'enzyme'},
                                    inputStyle={'margin-right':'5px'},
                                    className='enzyme-grid'
                                )
                            )
                        ], className='mb-3', style={'overflowY':'auto', 'maxHeight':'300px'}),
                        width=12
                    )
                ]),
            ]), 'CT', 'CodonTransformer')
    elif selected_model == 'JB':
        return (html.Div([
            html.Br(),
            html.Div([  
                html.Span("JohnBercow Options", className='model-options-title')
            ], style={'textAlign':'center'}),
            html.Br(),
            dbc.Row([
                dbc.Col([
                    html.H5('Order Name'),
                    dcc.Input(
                        id={'type':'param', 'name':'order_name'},
                        value='',
                        type='text',
                        placeholder='Your Order Name',
                        className="input-box",
                        style={'width':'100%'}
                    ),
                ], width=4),
                dbc.Col([
                    html.H5('Design Prefix'),
                    dcc.Input(
                        id={'type':'param', 'name':'design_prefix'},
                        type='text',
                        value='',
                        placeholder='Your design prefix',
                        className="input-box",
                        style={'width':'100%'}
                    )
                ], width=4),
                dbc.Col([
                    html.H5('Design ID'),
                    dcc.Input(
                        id={'type':'param', 'name':'design_id'},
                        placeholder='Your design ID',
                        value='',
                        className="input-box",
                        style={'width':'100%'}
                    )
                ], width=4)
            ]),
            html.Br(), #Separator
            dbc.Row([
                dbc.Col([
                    html.H5('Species'),
                    dcc.Dropdown(
                        id={'type':'param', 'name':'species_JB'},
                        options=[
                            {'label':'Escherichia coli', 'value':'e_coli'},
                            {'label':'Homo sapiens', 'value':'h_sapiens'}
                        ],
                        value='',
                    )
                ]),
                dbc.Col([
                    html.H5('Golden Gate Vector'),
                    dcc.Dropdown(
                        id={'type':'param', 'name':'golden_gate_vector_JB'},
                        options = get_gg_vectors(f'{script_dir}/utils/dna_extraction/JB/lab_vectors'),
                        value = ''
                    )
                ])
            ]),
            html.Br(), #Separator
            dbc.Row([
                dbc.Button('Other options', id='jb-options-button', n_clicks=0),
                dbc.Collapse(
                    html.Div(children=[
                        html.Br(), #Separator
                        dbc.Row([
                            dbc.Col([
                                html.H5('Skip IDT Optimization'),
                                dbc.RadioItems(
                                    id={'type':'param', 'name':'skip_idt_radio'},
                                    options=[
                                        {'label': 'True', 'value': True},
                                        {'label': 'False', 'value': False}
                                    ],
                                    value=False
                                )
                            ], width = 4),
                            dbc.Col([
                                html.H5('IDT Score Threshold'),
                                dbc.Input( #Add a callback that puts this in invalid state if the input is not a number
                                    id={'type':'param', 'name':'idt_score'},
                                    type='number',
                                    placeholder='IDT Score',
                                    value=7
                                )], width = 4                            
                            ),
                            dbc.Col([
                                html.H5('Sequence Max Length'),
                                dbc.Input(
                                    id={'type':'param', 'name':'sequence_max_length'},
                                    type='number',
                                    placeholder='Sequence Max Length',
                                    value=1500
                                )], width = 4
                            )
                        ]),
                        html.Br(), #Separator
                        dbc.Row([
                            dbc.Col([
                                html.H5('Starting Kmers Weight'),
                                dbc.Input(
                                    id={'type':'param', 'name':'starting_kmers_weight'},
                                    type='number',
                                    placeholder='Starting Kmers Weight',
                                    value=10
                                )], width=4
                            ),
                            dbc.Col([
                                html.H5('Number of Domesticator Steps'),
                                dbc.Input(
                                    id={'type':'param', 'name':'n_domesticator_steps'},
                                    type='number',
                                    placeholder='Number of Domesticator Steps',
                                    value=10
                                )], width=4
                            ),
                            dbc.Col([
                                 html.H5('Max Sequence Attempts'),
                                 dbc.Input(
                                     id={'type':'param', 'name':'max_attempts'},
                                     type='number',
                                     placeholder='Max Attempts',
                                     value=20
                                 )], width=4
                            ),
                        ]),
                        html.Br(), #Separator
                        dbc.Row([
                            dbc.Col([
                                html.H5('Print Heterooligomers'),
                                dbc.RadioItems(
                                    id={'type':'param', 'name':'print_heterooligomers'},
                                    options=[
                                        {'label': 'True', 'value': True},
                                        {'label': 'False', 'value': False}
                                    ],
                                    value=False
                                )]
                            ),
                            dbc.Col([
                                html.H5('No Layout'),
                                dbc.RadioItems(
                                    id={'type':'param', 'name':'no_layout'},
                                    options=[
                                        {'label': 'True', 'value': True},
                                        {'label': 'False', 'value': False}
                                    ],
                                    value=False
                                )]
                            ),
                            dbc.Col([
                                 html.H5('No Plasmids'),
                                 dbc.RadioItems(
                                     id={'type':'param', 'name':'no_plasmids'},
                                     options=[
                                         {'label': 'True', 'value': True},
                                         {'label': 'False', 'value': False}
                                     ],
                                     value=False
                                 )]
                            ),
                            dbc.Col([
                                html.H5('Verbose'),
                                dbc.RadioItems(
                                    id={'type':'param', 'name':'verbose'},
                                    options=[
                                        {'label': 'True', 'value': True},
                                        {'label': 'False', 'value': False}
                                    ],
                                    value=True
                                )]
                            ),
                            dbc.Col([
                                html.H5('Echo'),
                                dbc.RadioItems(
                                    id={'type':'param', 'name':'echo'},
                                    options=[
                                        {'label': 'True', 'value': True},
                                        {'label': 'False', 'value': False}
                                    ],
                                    value=False
                                )]
                            )
                        ])
                    ]), id='collapse-jb-options', is_open=False
                )
            ])
        ]), 'JB', 'JohnBercow')
    elif selected_model == 'PDB':
        return (html.Div(), 'PDB', 'PDB')

#Callback to open JB toggle
@callback(
    Output('collapse-jb-options', 'is_open'),
    Input('jb-options-button', 'n_clicks'),
    State('collapse-jb-options', 'is_open')
)
def toggle_filters_collapse_JB(n, is_open):
    if n:
        return not is_open
    return is_open

# Callback to extract hits
@callback(
    Output('extraction-viewer', 'children'),
    [
        Input('directory-dropdown', 'value'),
        Input('interval-component-extraction', 'n_intervals'),
        Input('execute-hits', 'n_clicks'),
        Input({'type':'param', 'name':dash.ALL}, 'value'),
        Input('extraction-selection', 'value')
    ],
    State('selected-model', 'data'),
    State('directory-dropdown', 'value'),
    prevent_initial_call=True
)

def extract_hits(working_dir, n, clicks, param_values, extraction_list, model_clicks, dir):
    import time
    #Create the hits folder
    hits_folder = os.path.join(working_dir, 'hits')
    fastas_folder = os.path.join(hits_folder, 'fastas')
    dnaseq_folder=os.path.join(hits_folder, 'dna_seqs')
    pdbs_folder = os.path.join(hits_folder, 'pdbs')
    if not os.path.exists(hits_folder):
        os.makedirs(hits_folder)
    if not os.path.exists(fastas_folder) or not os.path.exists(dnaseq_folder) or not os.path.exists(pdbs_folder):
        os.makedirs(fastas_folder, exist_ok=True)
        os.makedirs(dnaseq_folder, exist_ok=True)
        os.makedirs(pdbs_folder, exist_ok=True)

    dnaseq_folder = dnaseq_folder + '/'


    #Name of the multifastas output
    multifastas_output = 'all_sequences.fasta'
    multifastas_path = os.path.join(fastas_folder, multifastas_output)

    ctx = dash.callback_context

    # Message to the outfile viewer
    if ctx.triggered_id != "execute-hits":
        try:
            if model_clicks == "AF3":
                with open(f"{dir}/AF3_run/af3_logfile.log", 'r') as outfile:
                    stdout = outfile.read()
                    return stdout
            elif model_clicks == "PDB": 
                stdout = 'PDB OF THE HITS EXTRACTED. LIST OF PDBS EXTRACTED:\n' \
                     f'{os.listdir(pdbs_folder)}\n'
                return stdout
            elif model_clicks =="CT":
                return str("Waiting")
        except: 
            return str("waiting")

            #raise PreventUpdate

    else:
        #Getting the list of clicks
        model_states = ctx.states_list[0]

        chosen_model = model_states['value']

        print('CHOSEN MODEL:', chosen_model )

        # Getting all the parameters in a dict object, with name as key and value as value

        inputs_with_ids = ctx.inputs_list

        # Generator for dictionary
        params = (
            (str(param['id']['name']), param['value'])
            for item in inputs_with_ids
            if isinstance(item, list)
            for param in item
            if 'id' in param and 'name' in param['id'] and 'type' in param['id'] and 'param' in param['id']['type'] and 'value' in param
        )

        # Convert generator to dict when needed
        params_dict = dict(params)

        
        if chosen_model == 'CT':
            if extraction_list:
                for description in extraction_list:
            
                    #Get the fastas of the hits
        
                    extract_fasta_seq(os.path.join(pdbs_folder, description), fastas_folder)
            
            # Make the multifasta file

            multifastas(fastas_folder, multifastas_output)
            stdout = extract_dna_seq_CT(multifastas_path, dnaseq_folder, params_dict)
        
        elif chosen_model == "AF3":
            SCRIPT_DIR = os.path.dirname(__file__)
            subprocess.run(["mkdir -p "+str(dir)+"/AF3_run/slurm_logs"], shell=True)
            subprocess.run([f"python3 {SCRIPT_DIR}/BFmonitor/utils/AF3/AF3_runner.py --hits_dir {pdbs_folder} --directory {dir} > {dir}/AF3_run/af3_logfile.log"],shell=True)
            #AF3runner.main(pdbs_folder,dir)
            time.sleep(3)

            with open(f"{dir}/AF3_run/af3_logfile.log", 'r') as outfile:
                stdout = outfile.read()
        elif chosen_model == "JB":
            if extraction_list:
                for description in extraction_list:
                    # Get the fastas of the hits
                    extract_fasta_seq(os.path.join(pdbs_folder, description), fastas_folder)
                # Make the multifasta file
                multifastas(fastas_folder, multifastas_output)
                # Get the DNA sequences with JohnBercow
                stdout = extract_dna_seq_JB(multifastas_path, dnaseq_folder, params_dict)
                return stdout

        else:
            if extraction_list:
                for description in extraction_list:
                    # Get the PDBs
                    extract_pdbs(description,working_dir,pdbs_folder)
                    # Get the fasta sequence
                    extract_fasta_seq(os.path.join(pdbs_folder, description), fastas_folder)
                    # Add stats to the PDB file
                    add_stats_to_pdb(description, working_dir)
            stdout = 'PDB OF THE HITS EXTRACTED. LIST OF PDBS EXTRACTED:\n' \
                     f'{os.listdir(pdbs_folder)}\n' \

        # Make a csv file with the hits metrics
        
        metrics_df = pd.read_csv(os.path.join(working_dir, 'Scoring_Stats.csv'))
        filtered_df = metrics_df[metrics_df['description'].isin(extraction_list)]
        filtered_df.to_csv(f'{working_dir}/hits/hits_metrics.csv', index=False)
    

# Callback to toggle the Collapse for filters and axes in the sidebar
@callback(
    Output('filters-collapse', 'is_open'),
    Input('open-filters', 'n_clicks'),
    State('filters-collapse', 'is_open')
)

def toggle_filters_collapse(n, is_open):
    if n:
        return not is_open
    return is_open


# Callback for AF3 plot
 
@callback(
    Output("AF3-scatterplot", "figure"),
    Input("af3_x_axis", "value"),
    Input("af3_y_axis", "value"),
    State("directory-dropdown", "value")
)
def af3_plot(xaxis, yaxis, directory):
    import os
    af3_data_path = f'{directory}/Scoring_Stats_AF3.csv'
    if not os.path.exists(af3_data_path):
        return None # Not load figure, right now a placeholder
    scatterplot=scatter_plot_AF3(directory,xaxis,yaxis)
    return scatterplot

#Function to put the AF3 hits for processing
@callback(
    Output("af3-viewer", "children"),
    Input('interval-component-extraction', 'n_intervals'),
    State("directory-dropdown", "value")
    )
def af3_runner(intervals,dir):
    import os
    try:
        hits=""
        hits_path=os.path.join(str(dir)+"/hits/pdbs/","run_*")
        hits_files=glob.glob(hits_path, recursive=True)
        for f in hits_files:
            n=os.path.basename(f)
            hits=hits+str(">"+n+"\n")
    except:
        pass
    return hits

# end campaign
@callback(
    Input("stop-campaign","n_clicks"),
    prevent_initial_call=True
)
def stop_campaign(clk):
    #command = f'touch {working_dir}/campaign_done'
    #subprocess.run(command, shell=True)
    if not os.path.exists(f"{working_dir}/campaign_done"):
        #os.makedirs(f"{working_dir}/campaign_done")
        open(f"{working_dir}/campaign_done", 'a').close()


# Run the app
if __name__ == '__main__':
    app.run(debug=False, dev_tools_hot_reload = False, use_reloader=True,
                   host=hostname, port=port_number)


