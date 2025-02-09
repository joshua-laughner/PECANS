import nox

# To use uv to manage the venv, run as "nox --default-venv-backend uv"

@nox.session(python="3.13")
def test(session):
    session.install('.[all]')
    session.install('.')
    session.install('pytest')
    session.run('pytest')