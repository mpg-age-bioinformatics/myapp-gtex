from myapp import app, PAGE_PREFIX, PRIVATE_ROUTES
from flask_login import current_user
from flask_caching import Cache
import plotly.graph_objects as go
from flask import session, request
import dash
from dash import dcc, html
from dash.dependencies import Input, Output, State, MATCH, ALL
from dash.exceptions import PreventUpdate
from myapp.routes._utils import META_TAGS, navbar_A, protect_dashviews, make_navbar_logged
import dash_bootstrap_components as dbc
from myapp.routes.apps._utils import parse_import_json, parse_table, make_options, make_except_toast, ask_for_help, save_session, load_session, make_table, encode_session_app
import os
import sys
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
import fnmatch
from datetime import datetime
# import dash_table
import shutil
from time import sleep
from myapp import db
from myapp.models import UserLogging, PrivateRoutes
from ._app import read_menus, read_data, read_significant,  read_genes, gene_report, read_results_files, read_gene_expression, read_significant_genes, \
    filter_samples, filter_genes, filter_gene_expression, nFormat, read_dge,\
        make_volcano_plot, make_ma_plot, make_pca_plot, make_annotated_col, make_bar_plot, plot_height

from pyflaski.violinplot import make_figure, figure_defaults

        
PYFLASKI_VERSION=os.environ['PYFLASKI_VERSION']
PYFLASKI_VERSION=str(PYFLASKI_VERSION)


FONT_AWESOME = "https://use.fontawesome.com/releases/v5.7.2/css/all.css"

dashapp = dash.Dash("index",url_base_pathname=f'{PAGE_PREFIX}/', meta_tags=META_TAGS, server=app, external_stylesheets=[dbc.themes.BOOTSTRAP, FONT_AWESOME], title=app.config["APP_TITLE"], assets_folder=app.config["APP_ASSETS"])# , assets_folder="/flaski/flaski/static/dash/")

protect_dashviews(dashapp)

if app.config["SESSION_TYPE"] == "sqlalchemy":
    import sqlalchemy
    engine = sqlalchemy.create_engine(app.config["SQLALCHEMY_DATABASE_URI"] , echo=True)
    app.config["SESSION_SQLALCHEMY"] = engine
elif app.config["CACHE_TYPE"] == "RedisCache" :
    cache = Cache(dashapp.server, config={
        'CACHE_TYPE': 'RedisCache',
        'CACHE_REDIS_URL': 'redis://:%s@%s' %( os.environ.get('REDIS_PASSWORD'), app.config['REDIS_ADDRESS'] )  #'redis://localhost:6379'),
    })
elif app.config["CACHE_TYPE"] == "RedisSentinelCache" :
    cache = Cache(dashapp.server, config={
        'CACHE_TYPE': 'RedisSentinelCache',
        'CACHE_REDIS_SENTINELS': [ 
            [ os.environ.get('CACHE_REDIS_SENTINELS_address'), os.environ.get('CACHE_REDIS_SENTINELS_port') ]
        ],
        'CACHE_REDIS_SENTINEL_MASTER': os.environ.get('CACHE_REDIS_SENTINEL_MASTER')
    })

def change_table_minWidth(tb,minwidth):
    st=tb.style_table
    st["minWidth"]=minwidth
    tb.style_table=st
    return tb

def change_fig_minWidth(fig,minwidth):
    st=fig.style
    st["minWidth"]=minwidth
    fig.style=st
    return fig

def get_tables(cache,genders,tissues,groups,genenames,geneids):
    genes=read_genes(cache)
    data=read_data(cache)
    sigdf_=read_significant(cache)
    sigdf=sigdf_.drop(["file"],axis=1)

    if genders:
        data=data[data["gender"].isin(genders)]
        sigdf=sigdf[sigdf["gender"].isin(genders)]
        lgenders=len(genders)
    else:
        lgenders=0

    if tissues:
        data=data[data["tissue"].isin(tissues)]
        sigdf=sigdf[sigdf["tissue"].isin(tissues)]
        ltissues=len(tissues)
    else:
        ltissues=0

    if groups:
        data=data[ ( data["group_1"].isin(groups) ) | ( data["group_2"].isin(groups) )  ]
        sigdf=sigdf[ ( sigdf["group_1"].isin(groups) ) | ( sigdf["group_2"].isin(groups) )  ]

    if genenames or geneids :
        
        if genenames :
            lgenenames=len(genenames)
            sigdf_=sigdf[ ( sigdf["gene_name"].isin(genenames) ) ]
            genes_=genes[ ( genes["gene_name"].isin(genenames) ) ]
        else:
            lgenenames=0
            sigdf_=pd.DataFrame()
            genes_=pd.DataFrame()

        if geneids :
            lgeneids=len(geneids)
            sigdf__=sigdf[ ( sigdf["gene_id"].isin(geneids) ) ]
            genes__=genes[ ( genes["gene_id"].isin(geneids) ) ]
    
        else:
            lgeneids=0
            sigdf__=pd.DataFrame()
            genes__=pd.DataFrame()

        sigdf=pd.concat( [sigdf_, sigdf__ ] )
        sigdf=sigdf.drop_duplicates()

        genes=pd.concat( [genes_, genes__ ] )
        genes=genes.drop_duplicates()

    else:
        lgenenames=0
        lgeneids=0

    data=make_table(data,"data")
    sigdf=make_table(sigdf,"sigdf")

    if ( lgenders == 1 ) and ( ltissues == 1 ) and ( (lgenenames ==1 ) or (lgeneids == 1) ) :

        geneid=genes["gene_id"].tolist()[0]
        df=gene_report(cache, genders,tissues,geneid)
        df=df[["SAMPID","AGE","0","DTHHRDY", "SEX", "SMTS","SMTSD"]]
        df=df[2:]
        df["0"]=df["0"].astype(float)
        df=df.rename(columns={"0":"TPM"})
        df=df.sort_values(by=["AGE","SMTSD"],ascending=True)

        pa=figure_defaults()
        # session_file={"filename":"<from.gtex.app>", "last_modified":)}

        gene_name=genes["gene_name"].tolist()[0]
        gender=genders[0]
        tissue=tissues[0]

        pa["style"]="Violinplot and Swarmplot"
        pa['title']=f'{gene_name}, {tissue}, {gender}'
        pa["x_val"]="AGE"
        pa["y_val"]="TPM"
        pa["vals"]=[None]+df.columns.tolist()
        pa["xlabel"]="AGE"
        pa["ylabel"]="TPM"      

        session_data={ "session_data": {"app": { "violinplot": {"filename":"<from.gtex.app>" ,'last_modified':datetime.timestamp( datetime.now()),"df":df.to_json(),"pa":pa} } } }
        session_data["APP_VERSION"]=app.config['APP_VERSION']
        session_data["PYFLASKI_VERSION"]=PYFLASKI_VERSION

    else:

        df=None
        pa=None
        session_data=None

    return data, sigdf, df, pa, session_data

dashapp.layout=html.Div( 
    [ 
        dcc.Store( data=str(uuid.uuid4()), id='session-id' ),
        dcc.Location( id='url', refresh=False ),
        html.Div( id="protected-content" ),
    ] 
)

card_label_style={"margin-top":"5px"}
card_input_style={"width":"100%","height":"35px"}
card_body_style={ "padding":"2px", "padding-top":"4px"}


@dashapp.callback(
    Output('protected-content', 'children'),
    Input('session-id', 'data')
    )
def make_layout(session_id):
    if "gtex_" in PRIVATE_ROUTES :
        appdb=PrivateRoutes.query.filter_by(route="gtex_").first()
        if not appdb:
            return dcc.Location(pathname=f"{PAGE_PREFIX}/", id="index")
        allowed_users=appdb.users
        if not allowed_users:
            return dcc.Location(pathname=f"{PAGE_PREFIX}/", id="index")
        if current_user.id not in allowed_users :
            allowed_domains=appdb.users_domains
            if current_user.domain not in allowed_domains:
                return dcc.Location(pathname=f"{PAGE_PREFIX}/", id="index")
        # allowed_ips=app.config['WHITELISTED_IPS'].split(',') if app.config['WHITELISTED_IPS'] else []
        # user_ip=request.headers.get('X-Real-IP')
        # if allowed_ips and not any(fnmatch.fnmatch(user_ip, allowed_ip) for allowed_ip in allowed_ips):
        #     return dcc.Location(pathname=f"{PAGE_PREFIX}/", id="index")

    ## check if user is authorized
    eventlog = UserLogging(email=current_user.email, action="visit gtex")
    db.session.add(eventlog)
    db.session.commit()

    def make_loading(children,i):
        return dcc.Loading(
            id=f"menu-load-{i}",
            type="default",
            children=children,
        )

    protected_content=html.Div(
        [
            make_navbar_logged("GTEX",current_user),
            html.Div(id="app_access"),
            html.Div(id="redirect-violin"),
            dcc.Store(data=str(uuid.uuid4()), id='session-id'),
            dbc.Row(
                [
                    dbc.Col( 
                        [
                            dbc.Card(
                                [
                                    html.H5("Filters", style={"margin-top":10}), 
                                    html.Label('Gender'), make_loading( dcc.Dropdown( id='opt-genders', multi=True), 1), #opt-datasets
                                    html.Label('Tissue',style={"margin-top":10}),  make_loading( dcc.Dropdown( id='opt-tissues', multi=True), 3 ), #opt-samples
                                    html.Label('Age groups',style={"margin-top":10}),  make_loading( dcc.Dropdown( id='opt-groups', multi=True), 2 ), # opt-groups
                                    html.Label('Gene names',style={"margin-top":10}),  make_loading( dcc.Dropdown( id='opt-genenames', multi=True), 4 ),
                                    html.Label('Gene IDs',style={"margin-top":10}),  make_loading( dcc.Dropdown( id='opt-geneids', multi=True), 5 ),
                                    html.Label('Download file prefix',style={"margin-top":10}), 
                                    dcc.Input(id='download_name', value="gtex", type='text',style={"width":"100%", "height":"34px"})
                                ],
                                body=True
                            ),
                            dbc.Button(
                                'Submit',
                                id='submit-button-state', 
                                color="secondary",
                                n_clicks=0, 
                                style={"width":"100%","margin-top":"2px","margin-bottom":"2px"}#,"max-width":"375px","min-width":"375px"}
                            )
                        ],
                        sm=12,md=6,lg=4,xl=3,
                        align="top",
                        style={"padding":"0px","height": "100%",'overflow': 'scroll',"margin-bottom":"50px"} 
                    ),               
                    dbc.Col( 
                        dcc.Loading(
                            id="loading-output-2",
                            type="default",
                            children=[ html.Div(id="my-output")],
                            style={"margin-top":"50%","height": "100%"} 
                        ),
                        style={"height": "100%","width": "100%",'overflow': 'scroll',"margin-bottom":"50px"})
                ],
                align="start",
                justify="left",
                className="g-0",
                style={"height":"100%","width":"100%","overflow":"scroll"}
            ),
            navbar_A,
        ],
        style={"height":"100vh","verticalAlign":"center"}
    )
    return protected_content

#read_results_files
#read_genes


@dashapp.callback(
    Output(component_id='opt-genders', component_property='options'),
    Output(component_id='opt-tissues', component_property='options'),
    Output(component_id='opt-groups', component_property='options'),
    Output(component_id='opt-genenames', component_property='options'),
    Output(component_id='opt-geneids', component_property='options'),
    Input('session-id', 'data')
    )
def update_menus(session_id):
    menus=read_menus(cache)
    genders=make_options( menus["genders"] )
    tissues=make_options( menus["tissues"] )
    groups=make_options( menus["groups"] )

    genes=read_genes(cache)
    genenames=list(set(genes["gene_name"]))
    genenames=make_options(genenames)
    geneids=list(set(genes["gene_id"]))
    geneids=make_options(geneids)

    return genders, tissues, groups, genenames, geneids

@dashapp.callback(
    Output('my-output','children'),
    Input('session-id', 'data'),
    Input('submit-button-state', 'n_clicks'),
    State("opt-genders", "value"),
    State("opt-tissues", "value"),
    State("opt-groups", "value"),
    State("opt-genenames", "value"),
    State("opt-geneids", "value"),
    State('download_name','value'),
)
def update_output(session_id, n_clicks, genders, tissues, groups, genenames, geneids, download_name):

    # genes=read_genes(cache)
    # data=read_data(cache)
    # sigdf_=read_significant(cache)
    # sigdf=sigdf_.drop(["file"],axis=1)

    # if genders:
    #     data=data[data["gender"].isin(genders)]
    #     sigdf=sigdf[sigdf["gender"].isin(genders)]
    #     lgenders=len(genders)
    # else:
    #     lgenders=0

    # if tissues:
    #     data=data[data["tissue"].isin(tissues)]
    #     sigdf=sigdf[sigdf["tissue"].isin(tissues)]
    #     ltissues=len(tissues)
    # else:
    #     ltissues=0

    # if groups:
    #     data=data[ ( data["group_1"].isin(groups) ) | ( data["group_2"].isin(groups) )  ]
    #     sigdf=sigdf[ ( sigdf["group_1"].isin(groups) ) | ( sigdf["group_2"].isin(groups) )  ]

    # if genenames or geneids :
        
    #     if genenames :
    #         lgenenames=len(genenames)
    #         sigdf_=sigdf[ ( sigdf["gene_name"].isin(genenames) ) ]
    #         genes_=genes[ ( genes["gene_name"].isin(genenames) ) ]
    #     else:
    #         lgenenames=0
    #         sigdf_=pd.DataFrame()
    #         genes_=pd.DataFrame()

    #     if geneids :
    #         lgeneids=len(geneids)
    #         sigdf__=sigdf[ ( sigdf["gene_id"].isin(geneids) ) ]
    #         genes__=genes[ ( genes["gene_id"].isin(geneids) ) ]
    
    #     else:
    #         lgeneids=0
    #         sigdf__=pd.DataFrame()
    #         genes__=pd.DataFrame()

    #     sigdf=pd.concat( [sigdf_, sigdf__ ] )
    #     sigdf=sigdf.drop_duplicates()

    #     genes=pd.concat( [genes_, genes__ ] )
    #     genes=genes.drop_duplicates()

    # else:
    #     lgenenames=0
    #     lgeneids=0

    # data=make_table(data,"data")
    # sigdf=make_table(sigdf,"sigdf")

    swarmplot=[
        html.Div(
            [
                dbc.Row(
                    [
                        dbc.Col( 
                            [
                                html.Div(
                                    [ 
                                        f"Please select:",\
                                        html.Br(),\
                                        "1 gender",\
                                        html.Br(),\
                                        "1 tissue",\
                                        html.Br(),\
                                        "1 gene name or id",\
                                        html.Br()
                                    ],
                                    style={"textAlign":"center","overflow-wrap": "break-word"}
                                ),
                            ],
                            sm=9,md=7, lg=5, xl=5, 
                            align="top",
                            style={"textAlign":"center", "height": "100%", "margin-top":100  },
                        ),
                        navbar_A,
                    ],
                    justify="center",
                    style={"min-height": "100vh", "margin-bottom":"0px","margin-left":"5px","margin-right":"5px"}
                )
            ]
        )
    ]

    ## if only one gender, one tissue and one gene render a swarm plot of fpkms/tmps (as supplied by gtex)
    ## https://storage.googleapis.com/adult-gtex/bulk-gex/v8/rna-seq/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz

    # print(lgenders,ltissues,lgenenames,lgeneids)
    # print(( lgenders == 1 ) and ( ltissues == 1 ) and ( (lgenenames ==1 ) or (lgeneids == 1) ))

    data, sigdf, df, pa, session_data = get_tables(cache,genders,tissues,groups,genenames,geneids)

    # if ( lgenders == 1 ) and ( ltissues == 1 ) and ( (lgenenames ==1 ) or (lgeneids == 1) ):
    if pa:

        # geneid=genes["gene_id"].tolist()[0]
        # df=gene_report(cache, genders,tissues,geneid)
        # df=df[["SAMPID","AGE","0","DTHHRDY", "SEX", "SMTS","SMTSD"]]
        # df=df[2:]
        # df["0"]=df["0"].astype(float)
        # df=df.rename(columns={"0":"TPM"})
        # df=df.sort_values(by=["AGE","SMTSD"],ascending=True)

        # pa=figure_defaults()
        # session_file={"filename":"<from.gtex.app>", "last_modified":)}

        # gene_name=genes["gene_name"].tolist()[0]
        # gender=genders[0]
        # tissue=tissues[0]

        # pa["style"]="Violinplot and Swarmplot"
        # pa['title']=f'{gene_name}, {tissue}, {gender}'
        # pa["x_val"]="AGE"
        # pa["y_val"]="TPM"
        # pa["vals"]=[None]+df.columns.tolist()
        # pa["xlabel"]="AGE"
        # pa["ylabel"]="TPM"      

        # session_data={ "session_data": {"app": { "violinplot": {"filename":"<from.gtex.app>" ,'last_modified':datetime.timestamp( datetime.now()),"df":df.to_json(),"pa":pa} } } }
        # session_data["APP_VERSION"]=app.config['APP_VERSION']
        # session_data["PYFLASKI_VERSION"]=PYFLASKI_VERSION

        fig=make_figure(df,pa)
        fig_config={ 'modeBarButtonsToRemove':["toImage"], 'displaylogo': False}
        fig=dcc.Graph(figure=fig,config=fig_config,  id="graph")

        download_bar=html.Div( 
            [   
                dbc.Button(id='btn-download-values',className="me-1", n_clicks=0, children='Download values', style={"margin-top":4, 'background-color': "#5474d8", "color":"white"}),
                dcc.Download(id="download-values")
            ],
            style={'display': 'inline-block'}
        ) 

        send_to_violinplot=html.Div( 
            [
                dbc.Button(id='btn-violin-app', n_clicks=0, children='Violin plot', className="me-1",
                style={"margin-top":4, \
                    "margin-left":4,\
                    "margin-right":4,\
                    'background-color': "#5474d8", \
                    "color":"white",\
                    })
            ],
            style={'display': 'inline-block'}
        )

        message="Specific tissue information can be found on the original data. Please consider downloading the values and/or investigating the SMTSD column on the violin plot to understand the distribution of the different tissues."

        swarmplot=[fig, message, html.Div( [download_bar, send_to_violinplot]) ]


    minwidth=["Data","Significant"]
    minwidth=len(minwidth) * 150
    minwidth = str(minwidth) + "px"

    output=dcc.Tabs( 
        [ 
            dcc.Tab(
                dcc.Loading(
                    id="loading-output-2",
                    type="default",
                    children=[ data ],
                    style={"margin-top":"50%","height": "100%"} 
                ), 
                label="Data", id="tab-data",
                style={"margin-top":"0%"}
            ),
            dcc.Tab(
                dcc.Loading(
                    id="loading-output-3",
                    type="default",
                    children=[ sigdf ],
                    style={"margin-top":"50%","height": "100%"} 
                ), 
                label="Significant", id="tab-sig",
                style={"margin-top":"0%"}
            ),
            dcc.Tab(
                dcc.Loading(
                    id="loading-output-4",
                    type="default",
                    children=swarmplot,
                    style={"margin-top":"50%","height": "100%"} 
                ), 
                label="Swarmplot", id="tab-swarm",
                style={"margin-top":"0%"}
            )
        ]
    )

    return output


@dashapp.callback(
    Output("download-values", "data"),
    Input("btn-download-values", "n_clicks"),
    State("opt-genders", "value"),
    State("opt-tissues", "value"),
    State("opt-groups", "value"),
    State("opt-genenames", "value"),
    State("opt-geneids", "value"),
    State('download_name','value'),
    prevent_initial_call=True,
)
def download_values(n_clicks,genders, tissues, groups, genenames, geneids, download_name):

    data, sigdf, df, pa, session_data =get_tables(cache,genders,tissues,groups,genenames,geneids)

    # genes=read_genes(cache)

    # if genenames or geneids :
        
    #     if genenames :
    #         genes_=genes[ ( genes["gene_name"].isin(genenames) ) ]
    #     else:
    #         genes_=pd.DataFrame()

    #     if geneids :
    #         genes__=genes[ ( genes["gene_id"].isin(geneids) ) ]
    
    #     else:
    #         genes__=pd.DataFrame()

    #     genes=pd.concat( [genes_, genes__ ] )
    #     genes=genes.drop_duplicates()

    # geneid=genes["gene_id"].tolist()[0]
    # df=gene_report(cache, genders,tissues,geneid)
    # df=df[["SAMPID","AGE","0","DTHHRDY", "SEX", "SMTS","SMTSD"]]
    # df=df[2:]
    # df["0"]=df["0"].astype(float)
    # df=df.rename(columns={"0":"TPM"})
    # df=df.sort_values(by=["AGE","SMTSD"],ascending=True)

    fileprefix=secure_filename(str(download_name))
    filename="%s.xlsx" %fileprefix
    return dcc.send_data_frame(df.to_excel, filename, sheet_name="gtex", index=False)

@dashapp.callback(
    Output('redirect-violin','children'),
    Input('btn-violin-app', 'n_clicks'),
    State("opt-genders", "value"),
    State("opt-tissues", "value"),
    State("opt-groups", "value"),
    State("opt-genenames", "value"),
    State("opt-geneids", "value"),
)
def to_violin_app(n_clicks, genders, tissues, groups, genenames, geneids):
    if n_clicks:

        data, sigdf, df, pa, session_data =get_tables(cache,genders,tissues,groups,genenames,geneids)

        # genes=read_genes(cache)
        # data=read_data(cache)
        # sigdf_=read_significant(cache)
        # sigdf=sigdf_.drop(["file"],axis=1)

        # if genders:
        #     data=data[data["gender"].isin(genders)]
        #     sigdf=sigdf[sigdf["gender"].isin(genders)]
        #     lgenders=len(genders)
        # else:
        #     lgenders=0

        # if tissues:
        #     data=data[data["tissue"].isin(tissues)]
        #     sigdf=sigdf[sigdf["tissue"].isin(tissues)]
        #     ltissues=len(tissues)
        # else:
        #     ltissues=0

        # if groups:
        #     data=data[ ( data["group_1"].isin(groups) ) | ( data["group_2"].isin(groups) )  ]
        #     sigdf=sigdf[ ( sigdf["group_1"].isin(groups) ) | ( sigdf["group_2"].isin(groups) )  ]

        # if genenames or geneids :
            
        #     if genenames :
        #         lgenenames=len(genenames)
        #         sigdf_=sigdf[ ( sigdf["gene_name"].isin(genenames) ) ]
        #         genes_=genes[ ( genes["gene_name"].isin(genenames) ) ]
        #     else:
        #         lgenenames=0
        #         sigdf_=pd.DataFrame()
        #         genes_=pd.DataFrame()

        #     if geneids :
        #         lgeneids=len(geneids)
        #         sigdf__=sigdf[ ( sigdf["gene_id"].isin(geneids) ) ]
        #         genes__=genes[ ( genes["gene_id"].isin(geneids) ) ]
        
        #     else:
        #         lgeneids=0
        #         sigdf__=pd.DataFrame()
        #         genes__=pd.DataFrame()

        #     sigdf=pd.concat( [sigdf_, sigdf__ ] )
        #     sigdf=sigdf.drop_duplicates()

        #     genes=pd.concat( [genes_, genes__ ] )
        #     genes=genes.drop_duplicates()

        # else:
        #     lgenenames=0
        #     lgeneids=0

        # data=make_table(data,"data")
        # sigdf=make_table(sigdf,"sigdf")

        # geneid=genes["gene_id"].tolist()[0]
        # df=gene_report(cache, genders,tissues,geneid)
        # df=df[["SAMPID","AGE","0","DTHHRDY", "SEX", "SMTS","SMTSD"]]
        # df=df[2:]
        # df["0"]=df["0"].astype(float)
        # df=df.rename(columns={"0":"TPM"})
        # df=df.sort_values(by=["AGE","SMTSD"],ascending=True)

        # pa=figure_defaults()
        # session_file={"filename":"<from.gtex.app>", "last_modified":)}

        # gene_name=genes["gene_name"].tolist()[0]
        # gender=genders[0]
        # tissue=tissues[0]

        # pa["style"]="Violinplot and Swarmplot"
        # pa['title']=f'{gene_name}, {tissue}, {gender}'
        # pa["x_val"]="AGE"
        # pa["y_val"]="TPM"
        # pa["vals"]=[None]+df.columns.tolist()
        # pa["xlabel"]="AGE"
        # pa["ylabel"]="TPM"      

        # session_data={ "session_data": {"app": { "violinplot": {"filename":"<from.gtex.app>" ,'last_modified':datetime.timestamp( datetime.now()),"df":df.to_json(),"pa":pa} } } }
        # session_data["APP_VERSION"]=app.config['APP_VERSION']
        # session_data["PYFLASKI_VERSION"]=PYFLASKI_VERSION
                        
        session_data=encode_session_app(session_data)
        session["session_data"]=session_data

        from time import sleep
        sleep(2)

        return dcc.Location(pathname=f"{PAGE_PREFIX}/violinplot/", id="index")



@dashapp.callback(
        Output("navbar-collapse", "is_open"),
        [Input("navbar-toggler", "n_clicks")],
        [State("navbar-collapse", "is_open")])
def toggle_navbar_collapse(n, is_open):
    if n:
        return not is_open
    return is_open
