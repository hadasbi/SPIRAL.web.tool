#!/var/www/html/SPIRAL.web.tool/spiral_venv/bin/python3.10

import sys
import logging

sys.path.insert(0, '/var/www/html/SPIRAL.web.tool')
sys.path.insert(0, '/var/www/html/SPIRAL.web.tool/spiral_venv/python3.9/site-packages/')

# Set up logging
logging.basicConfig(stream=sys.stderr, level=logging.DEBUG)

activate_this = '/var/www/html/SPIRAL.web.tool/spiral_venv/bin/activate_this.py'
with open(activate_this) as file_:
    exec(file_.read(), dict(__file__=activate_this))

# Import and run the Flask app
from main import app as application