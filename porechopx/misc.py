class MyHelpFormatter(argparse.HelpFormatter):
    """
    This is a custom formatter class for argparse. It allows for some custom formatting,
    in particular for the help texts with multiple options (like bridging mode and verbosity level).
    http://stackoverflow.com/questions/3853722
    """
    def __init__(self, prog):
        terminal_width = shutil.get_terminal_size().columns
        os.environ['COLUMNS'] = str(terminal_width)
        max_help_position = min(max(24, terminal_width // 3), 40)
        super().__init__(prog, max_help_position=max_help_position)

    def _get_help_string(self, action):
        help_text = action.help
        if action.default != argparse.SUPPRESS and 'default' not in help_text.lower() and \
                action.default is not None:
            help_text += ' (default: ' + str(action.default) + ')'
        return help_text