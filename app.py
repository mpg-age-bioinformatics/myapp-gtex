from myapp import app, PAGE_PREFIX, PRIVATE_ROUTES
from flask_login import current_user
from flask_caching import Cache
import plotly.graph_objects as go
from flask import session
import dash
from dash import dcc, html
from dash.dependencies import Input, Output, State, MATCH, ALL
from dash.exceptions import PreventUpdate
from myapp.routes._utils import META_TAGS, navbar_A, protect_dashviews, make_navbar_logged
import dash_bootstrap_components as dbc
from myapp.routes.apps._utils import parse_import_json, parse_table, make_options, make_except_toast, ask_for_help, save_session, load_session, make_table, encode_session_app
import os
import uuid
# import traceback
import json
import base64
import pandas as pd
# import time
from werkzeug.utils import secure_filename
import humanize
from myapp.models import User
import stat
from datetime import datetime
# import dash_table
import shutil
from time import sleep
from myapp import db
from myapp.models import UserLogging, PrivateRoutes
from ._app import read_results_files, read_gene_expression, read_genes, read_significant_genes, \
    filter_samples, filter_genes, filter_gene_expression, nFormat, read_dge,\
        make_volcano_plot, make_ma_plot, make_pca_plot, make_annotated_col, make_bar_plot, plot_height
        

dashapp = dash.Dash("index",url_base_pathname=f'{PAGE_PREFIX}/', meta_tags=META_TAGS, server=app, external_stylesheets=[dbc.themes.BOOTSTRAP], title=app.config["APP_TITLE"], assets_folder=app.config["APP_ASSETS"])# , assets_folder="/flaski/flaski/static/dash/")

protect_dashviews(dashapp)

dashapp.layout=html.Div( [ 
                dcc.Location(id='url', refresh=False),
                html.Div(id="protected-content"),
                ] 
            )

@dashapp.callback(
    Output('protected-content', 'children'),
    Input('url', 'pathname'))
def make_layout(pathname):
    protected_content=html.Div(
        [
            make_navbar_logged("gtex",current_user),
            dbc.Container(
                dbc.Row( 
                    dbc.Col(
                        [
                            html.H1("Simple demo on how to add an app as an individual container.\nMy first changes.")
                        ],
                        align="center",
                    ),
                align="center",
                justify="center",
                style={'textAlign':'center',"height":"87vh"}
                ),
            ),
            navbar_A
        ],
        style={"height":"100vh","verticalAlign":"center"}
    )
    return protected_content


@dashapp.callback(
    Output("navbar-collapse", "is_open"),
    [Input("navbar-toggler", "n_clicks")],
    [State("navbar-collapse", "is_open")],
    )
def toggle_navbar_collapse(n, is_open):
    if n:
        return not is_open
    return is_open