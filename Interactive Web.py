
#Create Dash app
import dash
from dash import dcc, html
from dash.dependencies import Input, Output
import pandas as pd
import joblib
import numpy as np

# Load the model from the saved file
#model = joblib.load('knn_model_scaled_k5.pkl')

# Load clinical data for the correlation matrix (optional, not used in this section)
# clinical = pd.read_csv("Clinical data with mutation count.csv", sep=",")
#
# numeric_columns = clinical.select_dtypes(include=[np.number])
#
# # Compute the correlation matrix with numeric columns only
# correlation_matrix = numeric_columns.corr()
#
#
# correlation_with_target = numeric_columns.corr()['Neoplasm Disease Stage American Joint Committee on Cancer Code'].sort_values(ascending=False)
# print(correlation_with_target)
#
# #Print or log the correlation matrix (optional, for debugging)
# print(correlation_matrix)
#
# # correlation_matrix = clinical.corr()
# from dash import dcc, html, Input, Output
# import dash
# import numpy as np
# import pandas as pd
# import joblib
# from sklearn.neighbors import KNeighborsClassifier
#
# # Load the retrained model and scaler
# model = joblib.load('knn_model_retrained_subset.pkl')
# scaler = joblib.load('scaler.pkl')
#
# from dash import dcc, html, Input, Output
# import dash
# import numpy as np
# import pandas as pd
# import joblib
# from sklearn.neighbors import KNeighborsClassifier
#
# # Load the retrained model and scaler
# model = joblib.load('knn_model_retrained_subset.pkl')
# scaler = joblib.load('scaler.pkl')
#
# # Create the Dash app instance
# app = dash.Dash(__name__, suppress_callback_exceptions=True)
#
# # Layout for the app
# app.layout = html.Div(
#     children=[
#         html.H1("Thyroid Prediction Dataset", style={'color': 'black', 'text-align': 'center'}),
#         dcc.Tabs(
#             id='tabs',
#             value='tab-1',
#             children=[
#                 dcc.Tab(label='Survival Analysis Prediction', value='tab-1'),
#                 dcc.Tab(label='Data Visualization', value='tab-2'),
#                 dcc.Tab(label='Thyroid Cancer Stage Information', value='tab-3')
#             ],
#             colors={'border': 'pink', 'primary': 'pink', 'background': 'pink'}
#         ),
#         html.Div(id='tabs-content')
#     ],
#     style={'backgroundColor': 'white'}
# )
#
# @app.callback(
#     Output('tabs-content', 'children'),
#     [Input('tabs', 'value')]
# )
# def render_content(tab):
#     if tab == 'tab-1':
#         return html.Div([
#             html.H2("Survival Analysis Prediction"),
#             html.H3("Enter Patient Features for Prediction"),
#             html.Div([
#                 html.Label('Ethnicity Category:'),
#                 dcc.Input(id='ethnicity', type='text', placeholder='Enter Ethnicity (e.g., Hispanic, Non-Hispanic)'),
#                 html.Label('Diagnosis Age:'),
#                 dcc.Input(id='diagnosis_age', type='number', placeholder='Enter Age at Diagnosis'),
#                 html.Label('Race Category:'),
#                 dcc.Input(id='race', type='text', placeholder='Enter Race (e.g., White, Asian)'),
#                 html.Label('Sex:'),
#                 dcc.Input(id='sex', type='text', placeholder='Enter Sex (e.g., Male, Female)'),
#                 html.Label('Prior Cancer Diagnosis Occurrence:'),
#                 dcc.Input(id='Prior_Cancer', type='text', placeholder='Enter Prior Cancer Diagnosis (e.g., Yes or No)'),
#                 html.Label('Mutation Count:'),
#                 dcc.Input(id='Mutation', type='text', placeholder='Enter Mutation Count (e.g. Number)'),
#             ]),
#             html.Button('Submit', id='submit-button', n_clicks=0),
#             html.Hr(),
#             html.Div(id='prediction-output')
#         ])
#     elif tab == 'tab-2':
#         return html.Div([html.H2("Data Visualization"), html.P("Visualization will go here.")])
#     elif tab == 'tab-3':
#         return html.Div([
#             html.H2("Data Description"),
#             html.P(["This model uses clinical and genomic data from the ",
#                     html.A("GDC data portal", href="https://portal.gdc.cancer.gov/projects/TCGA-THCA", target="_blank"),
#                     " to predict thyroid cancer stages and survival analysis."]),
#             html.H3("Patient Information"),
#             html.P(["For more information on thyroid cancer staging and treatment, refer to the following resources:"]),
#             html.P([html.A("American Cancer Society",
#                            href="https://www.cancer.org/cancer/types/thyroid-cancer/detection-diagnosis-staging/staging.html",
#                            target="_blank")]),
#             html.P([html.A("Cancer.gov", href="https://www.cancer.gov/types/thyroid/patient/thyroid-treatment-pdq",
#                            target="_blank")])
#         ])
#     return html.Div('No content yet')
#
# @app.callback(
#     Output('prediction-output', 'children'),
#     [Input('submit-button', 'n_clicks')],
#     [Input('ethnicity', 'value'),
#      Input('diagnosis_age', 'value'),
#      Input('race', 'value'),
#      Input('sex', 'value'),
#      Input('Prior_Cancer', 'value'),
#      Input('Mutation', 'value')]
# )
# def make_prediction(n_clicks, ethnicity, diagnosis_age, race, sex, Prior_Cancer, Mutation):
#     if n_clicks > 0:
#         # Check if all inputs are provided
#         if None not in [ethnicity, diagnosis_age, race, sex, Prior_Cancer, Mutation]:
#             try:
#                 # Create input DataFrame for prediction
#                 input_df = pd.DataFrame([{
#                     'Ethnicity Category': ethnicity,
#                     'Diagnosis Age': diagnosis_age,
#                     'Race Category': race,
#                     'Sex': sex,
#                     'Prior Cancer Diagnosis Occurrence ': Prior_Cancer,
#                     'Mutation Count': Mutation
#                 }])
#
#                 # Encode categorical columns as necessary
#                 for col in input_df.select_dtypes(include=['object']).columns:
#                     input_df[col] = input_df[col].astype('category').cat.codes
#
#                 # Ensure columns match those used in training
#                 selected_features = ['Ethnicity Category', 'Diagnosis Age', 'Race Category', 'Sex',
#                                      'Prior Cancer Diagnosis Occurrence ', 'Mutation Count']
#                 input_df = input_df[selected_features]
#
#                 # Scale input data
#                 input_scaled = scaler.transform(input_df)
#
#                 # Make prediction
#                 prediction = model.predict(input_scaled)
#                 prediction_text = f"Predicted Survival Stage: {prediction[0]}"
#
#                 return html.Div([html.H4(prediction_text, style={'color': 'green'})])
#             except Exception as e:
#                 return html.Div([html.H4(f"Error during prediction: {str(e)}", style={'color': 'red'})])
#         else:
#             return html.Div([html.H4("Please fill in all the feature inputs.", style={'color': 'red'})])
#     return html.Div('')
#
# # Run the Dash app
# if __name__ == '__main__':
#     app.run_server(debug=True)

