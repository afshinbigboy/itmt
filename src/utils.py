import sys


class ColorPrint:

    @staticmethod
    def print_fail(*messages, end = '\n'):
        msg = ' '.join(str(m) for m in messages)
        sys.stderr.write('\x1b[1;31m' + msg.strip() + '\x1b[0m' + end)

    @staticmethod
    def print_pass(*messages, end = '\n'):
        msg = ' '.join(str(m) for m in messages)
        sys.stdout.write('\x1b[1;32m' + msg.strip() + '\x1b[0m' + end)

    @staticmethod
    def print_warn(*messages, end = '\n'):
        msg = ' '.join(str(m) for m in messages)
        sys.stderr.write('\x1b[1;33m' + msg.strip() + '\x1b[0m' + end)

    @staticmethod
    def print_info(*messages, end = '\n'):
        msg = ' '.join(str(m) for m in messages)
        sys.stdout.write('\x1b[1;34m' + msg.strip() + '\x1b[0m' + end)

    @staticmethod
    def print_bold(*messages, end = '\n'):
        msg = ' '.join(str(m) for m in messages)
        sys.stdout.write('\x1b[1;37m' + msg.strip() + '\x1b[0m' + end)