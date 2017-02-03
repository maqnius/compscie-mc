import tempfile
import os
import sys
import pytest
import shutil
import pkg_resources

test_pkg = 'particlesim'
cover_pkg = test_pkg

junit_xml = os.path.join(os.getenv('CIRCLE_TEST_REPORTS', '.'), 'junit.xml')

# chdir to an path outside of conda-bld, which is known to persist the build phase
run_dir = tempfile.mkdtemp()
os.chdir(run_dir)

pytest_args = ("-v --pyargs {test_pkg} "
               "--cov={cover_pkg} "
               "--cov-report=xml "
               "--junit-xml={junit_xml} "
               .format(test_pkg=test_pkg, cover_pkg=cover_pkg, junit_xml=junit_xml)
               .split(' '))
print("args:", pytest_args)
res = pytest.main(pytest_args)

# copy it to home, so we can process it with codecov etc.
shutil.copy('coverage.xml', os.path.expanduser('~/particlesim/'))

sys.exit(res)
