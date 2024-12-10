# Create Dash app
import dash
from dash import dcc, html
from dash.dependencies import Input, Output
import pandas as pd
import joblib
import numpy as np

# Load the model from the saved file
model = joblib.load('knn_model_retrained_subset.pkl')
scaler = joblib.load('scaler.pkl')

# Load clinical data for visualization
clinical = pd.read_csv("Clinical data with mutation count.csv", sep=",")

# Create the Dash app instance
app = dash.Dash(__name__, suppress_callback_exceptions=True)

# Layout for the app
app.layout = html.Div(
    children=[
        html.H1("Thyroid Prediction Dataset", style={'color': 'black', 'text-align': 'center'}),
        dcc.Tabs(
            id='tabs',
            value='tab-1',
            children=[
                dcc.Tab(label='Survival Analysis Prediction', value='tab-1'),
                dcc.Tab(label='Data Visualization', value='tab-2'),
                dcc.Tab(label='Thyroid Cancer Stage Information', value='tab-3')
            ],
            colors={'border': 'pink', 'primary': 'pink', 'background': 'pink'}
        ),
        html.Div(id='tabs-content')
    ],
    style={'backgroundColor': 'white'}
)

@app.callback(
    Output('tabs-content', 'children'),
    [Input('tabs', 'value')]
)
def render_content(tab):
    if tab == 'tab-1':
        return html.Div([
            html.H2("Survival Analysis Prediction"),
            html.H3("Enter Patient Features for Prediction"),
            html.Div([
                html.Label('Ethnicity Category:'),
                dcc.Input(id='ethnicity', type='text', placeholder='Enter Ethnicity (e.g., Hispanic, Non-Hispanic)'),
                html.Label('Diagnosis Age:'),
                dcc.Input(id='diagnosis_age', type='number', placeholder='Enter Age at Diagnosis'),
                html.Label('Race Category:'),
                dcc.Input(id='race', type='text', placeholder='Enter Race (e.g., White, Asian)'),
                html.Label('Sex:'),
                dcc.Input(id='sex', type='text', placeholder='Enter Sex (e.g., Male, Female)'),
                html.Label('Prior Cancer Diagnosis Occurrence:'),
                dcc.Input(id='Prior_Cancer', type='text', placeholder='Enter Prior Cancer Diagnosis (e.g., Yes or No)'),
                html.Label('Mutation Count:'),
                dcc.Input(id='Mutation', type='number', placeholder='Enter Mutation Count (e.g. Number)'),
            ]),
            html.Button('Submit', id='submit-button', n_clicks=0),
            html.Hr(),
            html.Div(id='prediction-output')
        ])
    elif tab == 'tab-2':
        # Dropdowns and plots for data visualization
        dataset_columns = clinical.columns.tolist()

        return html.Div([
            html.H2("Data Visualization"),
            html.Div([
                html.Label('Select X-axis:'),
                dcc.Dropdown(
                    id='column-dropdown_x',
                    options=[{'label': col, 'value': col} for col in dataset_columns],
                    value=dataset_columns[0],
                    style={'width': '48%', 'display': 'inline-block', 'margin-right': '2%'}
                ),
                html.Label('Select Y-axis:'),
                dcc.Dropdown(
                    id='column-dropdown_y',
                    options=[{'label': col, 'value': col} for col in dataset_columns],
                    value=dataset_columns[1],
                    style={'width': '48%', 'display': 'inline-block'}
                )
            ], style={'margin-bottom': '20px'}),
            dcc.Graph(id='scatter-plot'),
            dcc.Graph(id='histogram-plot')
        ])
    elif tab == 'tab-3':
        return html.Div([
            html.H2("Data Description"),
            html.P(["This model uses clinical and genomic data from the ",
                    html.A("GDC data portal", href="https://portal.gdc.cancer.gov/projects/TCGA-THCA", target="_blank"),
                    " to predict thyroid cancer stages and survival analysis."]),
            html.H3("Patient Information"),
            html.P(["For more information on thyroid cancer staging and treatment, refer to the following resources:"]),
            html.P([html.A("American Cancer Society",
                           href="https://www.cancer.org/cancer/types/thyroid-cancer/detection-diagnosis-staging/staging.html",
                           target="_blank")]),
            html.P([html.A("Cancer.gov", href="https://www.cancer.gov/types/thyroid/patient/thyroid-treatment-pdq",
                           target="_blank")])
        ])
    return html.Div('No content yet')

@app.callback(
    [Output('scatter-plot', 'figure'),
     Output('histogram-plot', 'figure')],
    [Input('column-dropdown_x', 'value'),
     Input('column-dropdown_y', 'value')]
)
def update_visualizations(x_column, y_column):
    if x_column and y_column:
        scatter_figure = {
            'data': [{
                'x': clinical[x_column],
                'y': clinical[y_column],
                'type': 'scatter',
                'mode': 'markers',
                'marker': {'color': 'blue'}
            }],
            'layout': {
                'title': f'Scatterplot: {x_column} vs {y_column}',
                'xaxis': {'title': x_column},
                'yaxis': {'title': y_column}
            }
        }

        histogram_figure = {
            'data': [{
                'x': clinical[x_column],
                'type': 'histogram',
                'marker': {'color': 'green'}
            }],
            'layout': {
                'title': f'Histogram of {x_column}',
                'xaxis': {'title': x_column},
                'yaxis': {'title': 'Count'}
            }
        }

        return scatter_figure, histogram_figure

    return {}, {}

@app.callback(
    Output('prediction-output', 'children'),
    [Input('submit-button', 'n_clicks')],
    [Input('ethnicity', 'value'),
     Input('diagnosis_age', 'value'),
     Input('race', 'value'),
     Input('sex', 'value'),
     Input('Prior_Cancer', 'value'),
     Input('Mutation', 'value')]
)
def make_prediction(n_clicks, ethnicity, diagnosis_age, race, sex, Prior_Cancer, Mutation):
    if n_clicks > 0:
        if None not in [ethnicity, diagnosis_age, race, sex, Prior_Cancer, Mutation]:
            try:
                input_df = pd.DataFrame([{
                    'Ethnicity Category': ethnicity,
                    'Diagnosis Age': diagnosis_age,
                    'Race Category': race,
                    'Sex': sex,
                    'Prior Cancer Diagnosis Occurrence ': Prior_Cancer,
                    'Mutation Count': Mutation
                }])

                for col in input_df.select_dtypes(include=['object']).columns:
                    input_df[col] = input_df[col].astype('category').cat.codes

                selected_features = ['Ethnicity Category', 'Diagnosis Age', 'Race Category', 'Sex',
                                     'Prior Cancer Diagnosis Occurrence ', 'Mutation Count']
                input_df = input_df[selected_features]

                input_scaled = scaler.transform(input_df)

                prediction = model.predict(input_scaled)
                prediction_text = f"Predicted Survival Stage: {prediction[0]}"

                return html.Div([html.H4(prediction_text, style={'color': 'green'})])
            except Exception as e:
                return html.Div([html.H4(f"Error during prediction: {str(e)}", style={'color': 'red'})])
        else:
            return html.Div([html.H4("Please fill in all the feature inputs.", style={'color': 'red'})])
    return html.Div('')

# Run the Dash app
if __name__ == '__main__':
    app.run_server(debug=True)
