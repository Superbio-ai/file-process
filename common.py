from pandas import DataFrame
from numpy import nan
from numpy import number as numpy_number


def create_tabular_response(data_df: DataFrame, include_feature_name=True):
    rows = data_df.round(2).replace({nan: None}).to_dict(orient='records')
    indices = list(data_df.index)
    if include_feature_name:
        for index, value in enumerate(indices):
            rows[index]['Feature Name'] = value
    return rows


def create_tabular_csv_response(data_df: DataFrame):
    numeric_columns = data_df.select_dtypes(include=numpy_number).columns
    rows = data_df.replace({nan: None})
    rows[numeric_columns] = rows[numeric_columns].round(2)
    rows = rows.to_dict(orient='records')
    return rows
