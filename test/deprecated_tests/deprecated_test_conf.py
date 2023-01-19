import pytest
from thermoengine import phases
from thermoengine import model

@pytest.fixture(scope="module")
def modelDB():
    return model.Database()

@pytest.fixture(scope="module")
def modelDB_Stix():
    return model.Database(database='Stixrude')

# @pytest.fixture(scope="module")
# def modelDB_HP():
#     return model.Database(database='HollandAndPowell')
