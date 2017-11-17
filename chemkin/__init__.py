import os

# Define path to XML files based on relative position to this module.
_PATH_XML_FILES = os.path.join(os.path.dirname(__file__), 'xml-files/')


def pckg_xml_path (file_name):
    """Creates full path to XML file located in xml-files directory."""
    if file_name[-4:] != '.xml':
        file_name += '.xml'
    return _PATH_XML_FILES + file_name
