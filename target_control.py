#!/usr/bin/env python

"""
Internal multiphase target control tool for internal target applications.

Mar 2015 by TwiNN

"""

import sys
from PyQt5.QtWidgets import QApplication
from velocitywindow import VelocityWindow


def main():
    app = QApplication(sys.argv)
    window = VelocityWindow()
    window.show()
    sys.exit(app.exec_())

if __name__ == '__main__':
    main()