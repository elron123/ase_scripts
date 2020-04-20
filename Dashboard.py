#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 20 12:36:43 2020

@author: elron
"""

# Import the packages needed
import pandas as pd
import numpy as np
import dash
from dash.dependencies import Input, Output
import dash_core_components as dcc
import dash_html_components as html
import plotly.graph_objs as go
import datetime
import json



# Read in the solvent data set
solv = pd.read_csv('angle_development.csv', index_col = 0)

# Defien a color dictionary
color = {"green": "rgb(149, 192, 66)", "blue": "rgb(0, 172, 214)", "grey": "rgb(156, 157, 159)"}



def scatter_plot(x_feature, y_feature, clickData = None):
    trace = go.Scatter(y=solv[y_feature], 
                      x=solv[x_feature], 
                      mode='markers', 
                      marker=dict(size=12,
                                  symbol = "circle",
                                  color = color["green"],
                                  line=dict(width=2,
                                            color="white",
                                           ),
                                 ),
                      text = solv[y_feature],
                      hovertemplate= "<b>%{text}</b><br>"+x_feature+": %{x}<br>"+y_feature+": %{y}<br>",
                      name = "",
                      selected = dict(marker = dict(color = color["blue"])),
                      unselected = dict(marker = dict(opacity = 0.5)),
                     )
    layout =go.Layout(title = dict(text = "Development of the Angle between various Atoms <br> for System [Zn] with Functional BP86", 
                                   x = 0.5,
                                  ),
                      xaxis = dict(title = "{}".format(x_feature), 
                                   #range = [solv[x_feature].min()*0.8, solv[x_feature].max()*1.2],
                                   linecolor = "black", 
                                   linewidth = 1,
                                   ticks = "outside",
                                   tickcolor = "black", 
                                   mirror = True,
                                   showgrid = True,
                                   gridcolor = color["grey"],
                                   gridwidth = 1,
                                   zeroline = True,
                                   zerolinecolor = color["grey"],
                                   zerolinewidth = 1,
                                  ),
                      yaxis = dict(title = "Angle [°]", 
                                   #range = [solv[y_feature].min()*0.8, solv[y_feature].max()*1.2],
                                   linecolor = "black", 
                                   linewidth = 1,
                                   ticks = "outside",
                                   tickcolor = "black", 
                                   mirror = True,
                                   showgrid = True,
                                   gridcolor = color["grey"],
                                   gridwidth = 1,
                                   zeroline = True,
                                   zerolinecolor = color["grey"],
                                   zerolinewidth = 1,
                                  ),    
                      font=dict(
                          family="Verdana, monospace",        
                          size=18,        
                          color="black"
                      ), 
                      paper_bgcolor = "white",
                      plot_bgcolor = "white",
                      hovermode = "closest",
                     )
    fig = go.Figure(data = trace, layout = layout)

    if clickData is not None: 
        fig.update_traces(
            overwrite=True,
            selectedpoints = [clickData["points"][0]["pointIndex"]]
        )
            
    return fig


#define the color scheme 
color = {"green": "rgb(149, 192, 66)", "blue": "rgb(0, 172, 214)", "grey": "rgb(156, 157, 159)"}

#define the border
border = None
#border = "thin grey solid"

app = dash.Dash()

#for each html element the style is defined
app.layout = html.Div([
    html.Div([
        html.H1(children = "[Zn]- with Functional BP86"),
    ], style = {"border": border, "text-align": "center", "margin-top": 20, "color": "black"}),
     html.Div([            
             html.Div([                
                     html.H2("Angle Development", 
                             style = {"text-align": "center", "margin-top": 58}),
             ], style = {"border": border, "margin": "1%",}),
        # Dropdown Menues x- and y- axis
             html.Div([
                 html.Div([                    
                     html.Label("Select a feature for the x-axis: "),                    
                     html.Div([                        
                         dcc.Dropdown(
                             id = "x_axis",
                             options  = [{"label": i, "value":i} for i in solv.columns],
                             value = "Step",),
                     ], style={"border": border, "margin-top": 10, "margin-bottom": 10}),
                 ], style = dict(border = border, width = "47%", display = "inline-block", margin = "1%")),
                 html.Div([ 
                     html.Label("Select a feature for the y-axis:"),
                     html.Div([                        
                         dcc.Dropdown(
                             id="y_axis",
                             options=[{"label": i, "value": i} for i in solv.columns],
                             value="Step",),
                     ], style={"margin-top": 10, "margin-bottom": 10}),
                 ], style = dict(border = border, width = "47%", float = "right", display = "inline-block", margin = "1%")),
             ], style = {"border": border, "margin": "1%"},
            ),
         dcc.Graph(id = "solvent_scatter", style = dict(border = border, margin = "1%"), 
                   figure={"layout":
                                    {"height":700
                                                            }
                                                                                                }),],),
#dropdown table
     
    ])

#The callback that updates the feature plot, when the dropdown is chosen or when a point is selected
@app.callback(
    dash.dependencies.Output("solvent_scatter", "figure"),
    [dash.dependencies.Input("x_axis", "value"),
     dash.dependencies.Input("y_axis", "value"),
     dash.dependencies.Input("solvent_scatter","clickData")]
    )
def update_scatter(x_axis,y_axis,clickData):
    plot = scatter_plot(x_axis, y_axis, clickData)
    return plot
     

       
app.run_server(port = 8052)