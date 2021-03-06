__author__ = 'Polychronis Patapis'
from PyQt4 import QtCore, QtGui
from core.pco_definitions import PixelFly
from threading import Thread
import os, time, sys, pickle
import pyqtgraph as pg
from astropy.io import fits
import numpy as np
import matlab.engine
from queue import Empty
import pygame, os, time, pickle
import win32api


class CameraWidget(QtGui.QWidget):
    """
    The CameraWidget class provides the user interface for the PCO PixelFly camera. It bases the connection to the
    camera through the pyPCOPixelFly.pco_definitions module. The basic framework of the class is PyQt4 an wrapper of the
    Qt framework and the pyqtgraph (url-here) module is essential for the use of this user interface.
     Dependencies:
     -- SC2_Cam.dll : the dynamic library that interfaces the camera hardware (please contain it in the same folder as
        the file).
     -- (Optional) App.ico : the application icon of pco (also needs to be in the same directory).

     Basic usage:
     Shortcuts:
     -- Ctrl + Q : Quits application
     -- Ctrl + R :Resets original scale to image
     Contact: Polychronis Patapis, patapisp@ethz.ch
    """

    def __init__(self, parent=None):
        QtGui.QWidget.__init__(self, parent)
        self.path = os.path.dirname(os.path.realpath("__file__"))
        self.save_dir = self.path
        self.camera = PixelFly(self.path)
        self.connected = False
        self.alive = False
        self.live_view_bool = False
        self.u = 1
        self.time_unit_dict = dict(us=1, ms=2)
        self.save_settings = self.load_settings()
        # set background color to dark gray
        self.setAutoFillBackground(True)
        p = self.palette()
        p.setColor(self.backgroundRole(), QtCore.Qt.darkGray)
        self.setPalette(p)

    def create_gui(self, MainWindow):
        """
        Creates user interface. Initializes all widgets of the application.
        :param MainWindow: The Main Application Window -> QtGui.MainWindow()
        :return:
        """
        # central widget of the Main Window
        self.central_widget = QtGui.QWidget(MainWindow)
        # set background color to dark gray
        self.central_widget.setAutoFillBackground(True)
        p = self.central_widget.palette()
        p.setColor(self.central_widget.backgroundRole(), QtCore.Qt.darkGray)
        self.central_widget.setPalette(p)
        # Grid layout to place all widgets
        self.widget_layout = QtGui.QGridLayout()
        # Graphics Layout Widget to put the image and histogram
        self.gw = pg.GraphicsLayoutWidget()
        # make margins around image items zero
        self.gw.ci.layout.setContentsMargins(0,0,0,0)
        # Graphics Layout Widget to put the crosscut curve plot
        self.gw_crosscut = pg.GraphicsLayoutWidget()

        MainWindow.setCentralWidget(self.central_widget)
        # the controls_layout contains all controls of the camera (eg. connection, exposure time, recording..)
        self.controls_layout = QtGui.QGridLayout()
        self.controls_layout.setSpacing(20)  # set spacing between widgets to 20 pixels
        # indicators_layout contains all indicators of the camera feed
        # The maximum count, the average count in the ROI region, buttons for ROI and crosscut, as well as
        # controls of the gray values if the image.
        self.indicators_layout = QtGui.QGridLayout()

        # ==============================================================================================================
        # CONTROL BUTTONS
        # ==============================================================================================================
        # Button to connect to the camera. Will turn red and display disconnect if it successfully connects.
        self.ConnectBtn = QtGui.QPushButton('CONNECT')
        self.controls_layout.addWidget(self.ConnectBtn, 0, 0)
        # layout for exposure time controls
        self.exsposure_time_layout = QtGui.QGridLayout()
        self.controls_layout.addItem(self.exsposure_time_layout, 2, 0, 4, 5)
        # 6 preset values of exposure time. They will be saved and reloaded through a python pickle file.
        preset_values = self.save_settings['exposure times']
        time_label1 = QtGui.QLabel("1")
        time_label2 = QtGui.QLabel("2")
        time_label3 = QtGui.QLabel("3")
        time_label4 = QtGui.QLabel("4")
        time_label5 = QtGui.QLabel("5")
        time_label6 = QtGui.QLabel("6")
        self.exp_time1 = QtGui.QPushButton(preset_values[0])
        self.exp_time2 = QtGui.QPushButton(preset_values[1])
        self.exp_time3 = QtGui.QPushButton(preset_values[2])
        self.exp_time4 = QtGui.QPushButton(preset_values[3])
        self.exp_time5 = QtGui.QPushButton(preset_values[4])
        self.exp_time6 = QtGui.QPushButton(preset_values[5])
        exposure_frame_title = QtGui.QLabel("Exposure time controls")
        self.exsposure_time_layout.addWidget(exposure_frame_title, 0, 0, 1, 3)
        self.exsposure_time_layout.addWidget(time_label1, 1, 0, 1, 1)
        self.exsposure_time_layout.addWidget(time_label2, 2, 0, 1, 1)
        self.exsposure_time_layout.addWidget(time_label3, 3, 0, 1, 1)
        self.exsposure_time_layout.addWidget(time_label4, 1, 2, 1, 1)
        self.exsposure_time_layout.addWidget(time_label5, 2, 2, 1, 1)
        self.exsposure_time_layout.addWidget(time_label6, 3, 2, 1, 1)
        self.exsposure_time_layout.addWidget(self.exp_time1, 1,1, 1, 1)
        self.exsposure_time_layout.addWidget(self.exp_time2, 2,1, 1, 1)
        self.exsposure_time_layout.addWidget(self.exp_time3, 3,1, 1, 1)
        self.exsposure_time_layout.addWidget(self.exp_time4, 1,3, 1, 1)
        self.exsposure_time_layout.addWidget(self.exp_time5, 2,3, 1, 1)
        self.exsposure_time_layout.addWidget(self.exp_time6, 3,3, 1, 1)
        # Edit line widget to input exposure time. It accepts us and ms units with the option of setting a float for
        # the ms time unit (eg. 1.5 ms)
        self.exp_time_in = QtGui.QLineEdit()
        # time units list
        self.time_units = QtGui.QComboBox()
        # save the time in one of the preset values.
        self.save_time = QtGui.QComboBox()

        self.exsposure_time_layout.addWidget(self.exp_time_in, 4, 2, 1, 3)
        self.exsposure_time_layout.addWidget(self.time_units, 4, 5, 1, 2)
        self.exsposure_time_layout.addWidget(self.save_time, 4, 0, 1, 2)

        # layout to host the recording controls
        self.recording_layout = QtGui.QGridLayout()
        self.controls_layout.addItem(self.recording_layout, 6, 0, 3, 3)
        recording_label = QtGui.QLabel("Recording controls")
        self.recording_layout.addWidget(recording_label, 0, 0, 1, 3)
        # Live button puts the camera in live view. Has to be stopped before exiting.
        self.LiveBtn = QtGui.QPushButton('LIVE')
        # Records the specified number of frames and lets the user name the file while adding 000x at the end
        # of the file name in FITS data format.
        self.RecordBtn = QtGui.QPushButton('RECORD')
        # stops live view/recording and disarms the camera
        self.StopBtn = QtGui.QPushButton('STOP')
        # Label for number of frames to save
        frame_lab = QtGui.QLabel('# frames to record:')
        # Edit line that accepts integers of the number of frames to save.
        self.FramesLab = QtGui.QLineEdit()
        self.recording_layout.addWidget(self.LiveBtn, 1, 0, 1, 1)
        self.recording_layout.addWidget(self.RecordBtn, 1, 1, 1, 1)
        #self.recording_layout.addWidget(self.StopBtn, 2, 0)
        self.recording_layout.addWidget(frame_lab, 2, 0, 1, 1)
        self.recording_layout.addWidget(self.FramesLab, 2, 1)

        # Callbacks for all the control buttons
        self.exp_time1.clicked.connect(self.exp_time_callback)
        self.exp_time2.clicked.connect(self.exp_time_callback)
        self.exp_time3.clicked.connect(self.exp_time_callback)
        self.exp_time4.clicked.connect(self.exp_time_callback)
        self.exp_time5.clicked.connect(self.exp_time_callback)
        self.exp_time6.released.connect(self.exp_time_callback)
        self.exp_time_list = [self.exp_time1, self.exp_time2, self.exp_time3, self.exp_time4,
                             self.exp_time5, self.exp_time6]
        # Add list options for time unit and save buttons.
        self.time_units.addItem("us")
        self.time_units.addItem("ms")
        self.time_units.activated[str].connect(self.onActivatedUnits)
        self.save_time.addItem("Save in")
        self.save_time.addItem("1")
        self.save_time.addItem("2")
        self.save_time.addItem("3")
        self.save_time.addItem("4")
        self.save_time.addItem("5")
        self.save_time.addItem("6")
        self.save_time.activated[str].connect(self.onActivatedSave)
        # Connect Enter/Return key press with callback for setting the exposure time.
        self.exp_time_in.returnPressed.connect(self.onReturnPress)
        # Connect callbacks for connect, live and stop buttons
        self.ConnectBtn.clicked.connect(self.connect_camera)
        self.ConnectBtn.setStyleSheet("background-color: darkCyan")

        self.FramesLab.setText('20')
        self.LiveBtn.clicked.connect(self.live_callback)
        #self.StopBtn.clicked.connect(self.stop_callback)
        self.RecordBtn.clicked.connect(self.record_callback)

        # layout to host the response matrix m
        self.ReponseMatrix_layout = QtGui.QGridLayout()
        self.controls_layout.addItem(self.ReponseMatrix_layout, 15, 0, 3, 3)
        # to start DM
        self.openDM_Btn = QtGui.QPushButton('open DM')
        # for response matrix
        self.ReponseMatrix_Btn = QtGui.QPushButton('RESPONSE MATRIX')
        # to close the DM after measurement
        self.closeDM_Btn = QtGui.QPushButton('close DM')
        # for the Zernike coefficient
        Enter_ZAmplitude = QtGui.QLabel("Amplitude:")
        self.Zernike_coef = QtGui.QLineEdit()
        self.Zernike_coef.setText("5")
        self.ReponseMatrix_layout.addWidget(Enter_ZAmplitude, 1, 0)
        self.ReponseMatrix_layout.addWidget(self.Zernike_coef, 1, 1)
        self.ReponseMatrix_layout.addWidget(self.ReponseMatrix_Btn, 3, 0)
        self.ReponseMatrix_layout.addWidget(self.closeDM_Btn, 2, 1)
        self.ReponseMatrix_layout.addWidget(self.openDM_Btn, 2, 0)
        self.ReponseMatrix_Btn.clicked.connect(self.Measure_ResponseMatrix)
        self.openDM_Btn.clicked.connect(self.open_DM)
        self.closeDM_Btn.clicked.connect(self.close_DM)

        # layout to host the SLM
        self.SLM_layout = QtGui.QGridLayout()
        self.controls_layout.addItem(self.SLM_layout, 10, 0, 3, 3)
        #to start DM
        self.activateSLM_Btn = QtGui.QPushButton('Initialize SLM')
        # to close the DM after measurement
       # self.closeSLM_Btn = QtGui.QPushButton('close SLM')
        # to create the phase map
        #self.createSLM_Btn = QtGui.QPushButton('create')
        # for the Zernike coefficient
        Enter_IoverD = QtGui.QLabel("file:")
        #self.Enter_IoverD .setText("D:\Xin LU\\real response matrix\\11X11\\Phase_shift_")
        self.file_address = QtGui.QLineEdit()
        self.file_address .setText("D:\Xin LU\\real response matrix\\9X9\\phase13_")
        self.SLM_layout.addWidget(Enter_IoverD, 1, 0)
        self.SLM_layout.addWidget(self.file_address , 1, 1)

        self.SLM_layout.addWidget(self.activateSLM_Btn, 2, 0)

        self.activateSLM_Btn.clicked.connect(self.activate_SLM)
        #self.closeSLM_Btn.clicked.connect(self.close_SLM)


        # ==============================================================================================================
        # IMAGE OPTIONS AND HANDLES
        # ==============================================================================================================
        # vb is a viewbox that contains the image item.
        self.vb = pg.ViewBox()
        # add the view box to the graphics layout
        self.gw.addItem(self.vb)
        # set the aspect while scaling to be locked, i.e. both axis scale the same.
        self.vb.setAspectLocked(lock=True, ratio=1)
        # invert Y axis -> PyQt <-> Numpy arrays convention
        self.vb.invertY()
        # Image Item is the image displaying item. Has a lot of options and the user can zoom in/out by pressing the
        # right mouse button and moving the mouse up/down. Furthermore by going over the image with the mouse will
        # indicate the coordinates and value.
        self.image = pg.ImageItem()
        self.vb.addItem(self.image)
        # Histogram of the displayed image. User can move the histogram axis and the gray values.
        self.hist = pg.HistogramLUTItem(self.image, fillHistogram=False)
        self.gw.addItem(self.hist)
        # initialize image container variable
        self.im = np.zeros((1392, 1040))
        # set image to display
        self.image.setImage(self.im)
        # set initial gray levels
        self.image.setLevels([200, 16383])
        self.hist.setHistogramRange(200, 16383)
        # Region Of Interest(ROI) widget that allows user to define a rectangle of tje image and the average count
        # within this will be displayed.
        #self.save_settings['ROI position']= ()
        self.roi = pg.ROI(pos=self.save_settings['ROI position'], size=self.save_settings['ROI size'])
        self.roi.addScaleHandle([1, 1], [0, 0])
        self.roi.alive = False
        self.vb.addItem(self.roi)
        self.roi.hide()
        # User can define line and place it on the image and the values profile will be plotted on the crosscut
        # graphics layout.
        self.line_roi = pg.LineSegmentROI([[680, 520], [720, 520]], pen='r')
        self.vb.addItem(self.line_roi)
        self.line_roi.hide()
        self.line_roi.alive = False
        # plot item to contain the crosscut curve
        crosscut_plot = pg.PlotItem()
        # crosscut curve that plot the data of the line
        self.crosscut_curve = pg.PlotCurveItem()
        self.gw_crosscut.addItem(crosscut_plot)
        crosscut_plot.addItem(self.crosscut_curve)
        self.gw_crosscut.hide()
        self.gw_crosscut.setFixedWidth(800)
        self.gw_crosscut.setFixedHeight(200)
        # make viewbox accept mouse hover events
        self.vb.acceptHoverEvents()
        # connect mouse moving event to callback
        self.vb.scene().sigMouseMoved.connect(self.mouseMoved)
        self.x, self.y = 0, 0  # mouse position
        # connect Ctrl + R key sequence to resetting the image to its original scale
        shortcut = QtGui.QShortcut(QtGui.QKeySequence('Ctrl+R'), MainWindow)
        shortcut.activated.connect(self.refresh_image)
        reset_btn = QtGui.QPushButton('Reset zoom')
        reset_btn.clicked.connect(self.refresh_image)
        # checkbox enabling log scale
        self.log_scale = QtGui.QCheckBox("Log scale")
        self.log_scale.stateChanged.connect(self.log_scale_callback)

        self.widget_layout.addWidget(self.gw, 0, 0, 6, 8)
        self.widget_layout.addWidget(self.gw_crosscut, 6, 3, 2, 6)
        self.widget_layout.addItem(self.controls_layout, 1, 8)
        self.widget_layout.addItem(self.indicators_layout, 7, 0, 2, 6)
        self.indicators_layout.addWidget(reset_btn, 2, 6, 1, 1)
        self.indicators_layout.addWidget(self.log_scale, 2, 7, 1, 1)
        # Indicator showing maxvalue of image being displayed
        self.max_indicator_lab = QtGui.QLabel('Max value')
        font = QtGui.QFont("Calibri", 18)
        self.max_indicator_lab.setFont(font)
        self.indicators_layout.addWidget(self.max_indicator_lab, 0,0,1,1)
        self.max_indicator = QtGui.QLabel(str(np.max(self.im)))
        self.max_indicator.setFont(font)
        self.indicators_layout.addWidget(self.max_indicator, 0,1,1,1)

        # Indicator showing average value within roi if it's selected
        self.roi_indicator = QtGui.QLabel('-')
        self.roi_indicator.setFont(QtGui.QFont("Calibri", 18))
        roi_indicator_lab = QtGui.QLabel('ROI average counts:')
        roi_indicator_lab.setFont(QtGui.QFont("Calibri", 18))
        self.indicators_layout.addWidget(roi_indicator_lab, 1, 0, 1, 1)
        self.indicators_layout.addWidget(self.roi_indicator, 1, 1, 1, 1)
        # Edit widget that allow setting the gray-levels
        self.gray_max = 16383
        self.gray_min = 200
        self.gray_max_edit = QtGui.QLineEdit(str(self.gray_max))
        self.gray_min_edit = QtGui.QLineEdit(str(self.gray_min))
        self.gray_min_lab = QtGui.QLabel('Min:')
        self.gray_max_lab = QtGui.QLabel('Max:')
        self.gray_min_edit.returnPressed.connect(self.set_gray_min)
        self.gray_max_edit.returnPressed.connect(self.set_gray_max)

        self.indicators_layout.addWidget(self.gray_min_lab, 2, 2, 1, 1)
        self.indicators_layout.addWidget(self.gray_max_lab, 2, 4, 1, 1)
        self.indicators_layout.addWidget(self.gray_min_edit, 2, 3, 1, 1)
        self.indicators_layout.addWidget(self.gray_max_edit, 2, 5, 1, 1)

        # Buttons for ROI and crosscut line
        roi_button = QtGui.QPushButton('ROI')
        crosscut_button = QtGui.QPushButton('Crosscut')
        self.indicators_layout.addWidget(roi_button, 2, 0, 1, 1)
        self.indicators_layout.addWidget(crosscut_button, 2, 1, 1, 1)
        roi_button.clicked.connect(self.roi_clicked)
        crosscut_button.clicked.connect(self.crosscut_clicked)
        #########################################
        self.central_widget.setLayout(self.widget_layout)
        # ==============================================================================================================
        # MENU BAR
        # ==============================================================================================================
        self.menubar = QtGui.QMenuBar(MainWindow)
        #self.menubar.setGeometry(QtCore.QRect(0, 0, 1027, 35))
        filemenu = self.menubar.addMenu('&File')

        exitAction = QtGui.QAction(QtGui.QIcon('exit.png'), '&Exit', self)
        exitAction.setShortcut('Ctrl+Q')
        exitAction.triggered.connect(self.closeEvent)
        filemenu.addAction(exitAction)
        MainWindow.setMenuBar(self.menubar)
        # ==============================================================================================================
        # STATUS BAR
        # ==============================================================================================================
        self.statusbar = QtGui.QStatusBar(MainWindow)
        font2 = QtGui.QFont("Calibri", 15)
        #self.statusbar.setGeometry(QtCore.QRect(0, 600, 1027, 35))
        self.statusbar.setStyleSheet("background-color: darkCyan")
        self.connection_status_lab = QtGui.QLabel('Connection status: ')
        self.connection_status_lab.setFont(font2)
        self.connection_status = QtGui.QLabel('Disconnected ')
        self.connection_status.setFont(font2)
        self.statusbar.addPermanentWidget(self.connection_status_lab)
        self.statusbar.addPermanentWidget(self.connection_status)
        self.display_status_lab = QtGui.QLabel('Display status: ')
        self.display_status_lab.setFont(font2)
        self.display_status = QtGui.QLabel('Idle ')
        self.display_status.setFont(font2)
        self.statusbar.addPermanentWidget(self.display_status_lab)
        self.statusbar.addPermanentWidget(self.display_status)
        self.measurement_status_lab = QtGui.QLabel('Measurement status: ')
        self.measurement_status_lab.setFont(font2)
        self.measurement_status = QtGui.QLabel(' - ')
        self.measurement_status.setFont(font2)
        self.statusbar.addPermanentWidget(self.measurement_status_lab)
        self.statusbar.addPermanentWidget(self.measurement_status)

        self.mouse_pos_lab = QtGui.QLabel('Mouse position: ')
        self.mouse_pos_lab.setFont(font2)
        self.mouse_pos = QtGui.QLabel(' - ')
        self.mouse_pos.setFont(font2)
        self.statusbar.addPermanentWidget(self.mouse_pos_lab)
        self.statusbar.addPermanentWidget(self.mouse_pos)
        MainWindow.setStatusBar(self.statusbar)

    def log_scale_callback(self):
        if self.log_scale.isChecked():
            self.image.setLevels([np.log(200), np.log(16383)])
            self.hist.setHistogramRange(np.log(200), np.log(16383))
        else:
            self.image.setLevels([200, 16383])
            self.hist.setHistogramRange(200, 16383)
        return

    def refresh_image(self):
        """
        Shortcut callback. If Ctrl+R is pressed the image scales back to its original range
        :return:
        """
        self.vb.autoRange()
        self.image.update()
        return

    def load_settings(self):
        """
        Load settings from previous session stored in gui_settings.p
        :return:
        """
        fname = self.path + '\\pco_settings.p'
        if os.path.isfile(fname):
            return pickle.load(open(fname, 'rb'))
        else:
            sets = {'ROI position': [696, 520], 'ROI size': 50, 'line position': [[10, 64], [120, 64]],
                    'exposure times': ['500 us', '800 us', '1 ms', '10 ms', '50 ms', '100 ms']}
            return sets

    def save_settings_return(self):
        """
        Save settings before exiting application
        :return:
        """
        fname = self.path + '\\pco_settings.p'
        times = []
        for btn in self.exp_time_list:
            times.append(btn.text())
        self.save_settings['exposure times'] = times
        self.save_settings['ROI position'] = self.roi.pos()
        self.save_settings['ROI size'] = self.roi.size()
        pickle.dump(self.save_settings, open( fname, "wb" ) )
        return


    def roi_clicked(self):
        """
        Callback to press of the ROI button. A rectangular roi will appear on the image corner.
        If active the roi will disappear.
        :return:
        """
        if self.roi.alive:
            self.roi.alive = False
            self.roi_indicator.setText('-')
            self.roi.hide()
        else:
            self.roi.alive = True
            self.roi.show()
        return

    def crosscut_clicked(self):
        """
        Callback to press of the line crosscut button. A line roi will appear on the image corner.
        If active the roi will disappear. The crosscut curve will also appear.
        :return:
        """
        if self.line_roi.alive:
            self.line_roi.alive = False
            self.gw_crosscut.hide()
            self.line_roi.hide()
        else:
            self.line_roi.alive = True
            self.gw_crosscut.show()
            self.line_roi.show()
        return

    def mouseMoved(self, event):
        """
        Mouse move callback. It displays the position and value of the mouse on the image on the statusbar, in
        the right corner.
        :param event: Mouse move event
        :return:
        """
        point = self.image.mapFromScene(event)
        self.x = int(point.x())
        self.y = int(point.y())
        # return if position out of image bounds
        if self.x < 0 or self.y < 0 or self.x > 1392 or self.y > 1040:
            return
        try:
            val = int(self.im[self.x, self.y])
            self.mouse_pos.setText('%i , %i : %i'%(self.x, self.y, val))
        except:
            pass
        return

    def roi_value(self):
        """
        Get data from ROI region and calculate average. The value will be displayed in the
        roi indicator label.
        :return:
        """
        data = self.roi.getArrayRegion(self.im, self.image)
        data_a, data_m, data_std = int(np.average(data)), int(np.max(data)), int(np.std(data))
        self.roi_indicator.setText('%i, Max: %i, STD: %i'%(data_a, data_m, data_std))
        return

    def line_roi_value(self):
        """
        Get data from line crosscut and plot them in the crosscut curve.
        :return:
        """
        data = self.line_roi.getArrayRegion(self.im, self.image)
        x_data = np.array(range(len(data)))
        self.crosscut_curve.setData(x_data, data)
        return

    def set_gray_max(self):
        """
        Set max value of graylevel. For the 14bit image the value is held up to 16383 counts.
        :return:
        """
        val = self.gray_max_edit.text()
        try:
            self.gray_max = int(val)
            if self.gray_max > 16383:
                self.gray_max = 16383
                self.gray_max_edit.setText('16383')
            self.image.setLevels([self.gray_min, self.gray_max])
            self.hist.setHistogramRange(self.gray_min, self.gray_max)
        except ValueError:
            pass
        return

    def set_gray_min(self):
        """
        Set min value of graylevel. For the 14bit image the value is held down to 0 counts.
        :return:
        """
        val = self.gray_min_edit.text()
        try:
            self.gray_min = int(val)
            if self.gray_min < 0:
                self.gray_min = 0
                self.gray_min_edit.setText('0')
            self.image.setLevels([self.gray_min, self.gray_max])
            self.hist.setHistogramRange(self.gray_min, self.gray_max)
        except ValueError:
            pass
        return

    def closeEvent(self, event):
        """
        Callback when exiting application. Ensures that camera is disconnected smoothly.
        :return:
        """
        if self.live_view_bool or self.alive:
            self.stop_callback()
        if self.connected:
            self.connect_camera()
        self.save_settings_return()
        QtGui.QApplication.closeAllWindows()
        QtGui.QApplication.instance().quit()
        
        return

    def onActivatedUnits(self, text):
        self.u = self.time_unit_dict[text]
        return

    def onActivatedSave(self, text):
        if text == "Save in":
            return
        which = int(text[-1])-1
        what = str(self.t) + ' ' + self.time_units.currentText()
        self.exp_time_list[which].setText(what)
        return

    def onReturnPress(self):
        text = self.exp_time_in.text()
        t, u = 0, 0
        try:
            if '.' in text or ',' in text and self.u == 2:
                t = int(float(text)*1000)
                u = 1
                self.t = float(text)
            else:
                self.t = int(text)
                t = self.t
                u = self.u             
            self.camera.exposure_time(t, u)
        except ValueError:
            pass
        return

    def connect_camera(self):
        """
        Connect to camera. If camera connection returns error report it and
        set connected status to False.
        :return:
        """
        if self.connected:
            err = self.camera.close_camera()
            self.connected = False
            self.ConnectBtn.setText('CONNECT')
            self.ConnectBtn.setStyleSheet("background-color: darkCyan")
            self.connection_status.setText('Disconnected')
        else:
            err = self.camera.open_camera()
            if not err:
                self.connection_status.setText('Error with connection')
                return
            self.connected = True
            self.ConnectBtn.setText('DISCONNECT')
            self.ConnectBtn.setStyleSheet("background-color: green")
            self.connection_status.setText('Connected')
            try:
                t, u = self.camera.get_exposure_time()                
                self.exp_time_in.setText(str(t))
                index = self.time_units.findText(u)
                if index >= 0:
                    self.u = self.time_unit_dict[u]
                    self.time_units.setCurrentIndex(index)
            except:
                pass
        return

    def exp_time_callback(self):
        """
        Set exposure time
        :param event: button press event
        :return:
        """
        which = self.sender().text()
        t, unit = which.split(sep=' ')
        unit_initial = unit
        try:
            if ('.' in t) or (',' in t) and (unit == 'ms'):
                self.t = int(float(t)*1000)
                unit = 'us'
            else:
                self.t = int(t)

            unit = self.time_unit_dict[unit]
            self.camera.exposure_time(self.t, unit)
            self.u = self.time_unit_dict[unit_initial]
            self.exp_time_in.setText(str(t))
            index = self.time_units.findText(unit_initial)
            if index >= 0:
                self.time_units.setCurrentIndex(index)
        except:
            pass
        return

    def live_callback(self):
        """
        Starts live view thread
        :return:
        """
        if self.connected:
            if self.alive:
                self.stop_callback()
                self.LiveBtn.setStyleSheet('background-color: lightGray')
                self.LiveBtn.setChecked(False)
            else:
                self.alive = True
                self.LiveBtn.setStyleSheet('background-color: darkCyan')
                self.live_thread = Thread(target=self.live_thread_callback)
                self.live_thread.setDaemon(True)
                self.live_thread.start()
                QtCore.QTimer.singleShot(500, self.update_image)
        else:
            self.display_status.setText('Error with live display')
        return

    def live_thread_callback(self):
        """
        Callback for thread that read images from buffer
        :return:
        """
        try:
            # Arm camera
            self.camera.arm_camera()
            print('Camera armed')
            # Set recording status to 1
            self.camera.start_recording()
            # Allocate buffers, default=2 buffers
            self.camera.allocate_buffer(3)
            self.camera._prepare_to_record_to_memory()
            self.display_status.setText('Live view.')
            self.record_live_thread = Thread(target=self.camera.record_to_memory_2)
            print('record thread created')
            self.record_live_thread.setDaemon(True)
            self.record_live_thread.start()
            print('record thread started')
            self.live_view_bool = True
            """
            Remember to manage all exceptions here. Look it up in PixelFly() class
            """
        except:
            self.stop_callback()
        return

    def open_DM(self):
        """
        open DM for the response matrix measurement
        if self.driver_info['UBS_pointer']<0 , keeping doing so until it is >0
        :return:
        """
        self.eng = matlab.engine.start_matlab()
        # start the DM
        (self.error_code, self.driver_info) = self.eng.openDM(6.0, nargout=2)
        print(self.driver_info)

        if self.error_code != 0:
            print('DM error')
            self.eng.quit()
        elif self.driver_info['USB_pointer'] <= 0:
            print('USB_pointer is invalid!!')
            self.eng.quit()

    def activate_SLM(self):
        # get path that file is saved in
        self.path = os.path.dirname(os.path.realpath("__file__"))
        # find connected screens
        screens = win32api.EnumDisplayMonitors()
        if len(screens) < 2:
            raise UserWarning('No second screen connected')
        self.hdc = screens[1][0]  # handle of second screen PYHANDLE
        self.dim = screens[1][2]  # screen dimensions
        self.left = self.dim[0]
        self.top = self.dim[1]
        self.right = self.dim[2]
        self.bottom = self.dim[3]
        self.width = abs(self.right - self.left)
        self.height = abs(self.bottom - self.top)
        if self.width == 1024:
            self.SLM_type = "LC-2012"
        else:
            self.SLM_type = "Pluto"
        self.size = (self.width, self.height)
        self.dimensions = (self.width, self.height, 3)
        # set Windows ENVIRONMENTAL VARIABLE of SDL window position to the top-left corner
        # of the second screen
        os.environ['SDL_VIDEO_WINDOW_POS'] = "%d,%d" % (self.left, self.top)
        pygame.display.init()
        self.screen = pygame.display.set_mode(self.size)
        # create surface object
        self.surface = pygame.surfarray.make_surface(np.zeros(self.dimensions, dtype=np.uint8))
        self.screen.blit(self.surface, (0, 0))
        pygame.display.flip()
        # read generated phase map
        self.SLM_map = []
        for index in range(4):
            file_name = self.file_address.text()+str(index)+'.fits'
            self.SLM_map.append(self.import_phase_map(file_name))

    def import_phase_map(self, file):
        """
        Imports an phase map stored in a file and saves it in an
        numpy array for use on the SLM
        :param file: A .txt or .fits file of a phase map
        :return: write an 8bit array of the correct dimensions in maps dict
        """
        hdu = fits.open(file)
        p = hdu[0].data
        m = np.zeros((self.width, self.height, 3), dtype=np.uint8)
        m[:, :, 0] = p
        m[:, :, 1] = p
        m[:, :, 2] = p
        return m

    def draw(self, p):
        """
        Draw array onto second screen.
        :param p: p is an uint8 numpy array with the correct dimensions
        :return:
        """
        surface = pygame.surfarray.make_surface(p)
        self.screen.blit(surface, (0, 0))
        pygame.display.flip()
        return

    def update_image(self):
        """
        Takes images from camera queue and displays them. If roi or crosscut is enabled, it updates the
        respective values/plot. The consumer loop works using the QtCore.QTimer.singleShot() method,
        that fires the function every x ms until it is interrupted.
        :return:
        """
        if not self.alive and (not self.camera.armed):
            if self.live_view_bool:
                self.live_view_bool = False
                time.sleep(0.1)
                self.record_live_thread.join()
                self.live_thread.join()
                del self.record_live_thread
                del self.live_thread
            self.display_status.setText('Idle')
            return
        try:
            # get newest frame from queue. Transpose it so that is fits the coordinates convention
            im = self.camera.q.get().T
            # check for log scale display
            if self.log_scale.isChecked():
                self.im = np.log(im)
            else:
                self.im = im
            try:
                # get max value from queue
                max_val = self.camera.q_m.get()
                self.max_indicator.setText(str(max_val))
            except Empty:
                pass
            # set new image data, with options autoLevels=False so that it doesn't change the grayvalues
            # autoRange=False so that it stays at the zoom level we want and autoHistogram=False so that it does
            # not interfere with the axis of the histogram
            self.image.setImage(self.im, autoLevels=False, autoRange=False, autoHistogramRange=False,
                                autoDownsample=True)
            # if roi button is clicked
            if self.roi.alive:
                self.roi_value()
            # if crosscut line is clicked
            if self.line_roi.alive:
                self.line_roi_value()
            # mouse position value update

            if 0 <= self.x <= 1392 and 0 <= self.y <= 1040:
                val = im[self.x, self.y]
                self.mouse_pos.setText('%i , %i : %.1f' % (self.x, self.y, val))

            # update image. Don't know if this is necessary..
            self.image.update()
        except Empty:
            pass
        # Run single shot timer again
        QtCore.QTimer.singleShot(20, self.update_image)


    def Measure_ResponseMatrix(self):
        """
        Record the images of the first 15 Zernike terms with Zernike_coefficient = self.Zernike_coef
        :return:
        """
        if not self.connected:
            return
        # check if camera had already been armed
        if self.alive and self.live_view_bool:
            try:
                self.live_view_bool = False
                self.camera.live = False  # stop loop that is producing frames
                #time.sleep(1)
                self.record_live_thread.join()
                self.live_thread.join()
                del self.record_live_thread
                del self.live_thread
            except:
                pass
        elif not self.alive:
            self.alive = True
            self.camera.arm_camera()
            self.camera.start_recording()
            self.camera.allocate_buffer()
            self.camera._prepare_to_record_to_memory()
            self.alive = True
        else:
            pass

        try:
            num_of_frames = int(self.FramesLab.text())
            self.Z_coef = int(self.Zernike_coef.text())
        except ValueError:
            return


        # probably there is a conflict between the camera and DM USB connections, so the 'USB_pointer' in the self.driver_info
        # is negative, we need to read the self.driver_info from DM_GUI
        #self.driver_info = np.load('current_DMdriver_info.npy').item()

        RM_filename = QtGui.QFileDialog.getSaveFileName(self, 'Save as..', self.save_dir)
        self.save_dir = os.path.dirname(RM_filename)

        for Zernike_term in range(1,28):   # 11 for first 10 Zernike terms

            #Zernike_coefficient = [0] * 15 # initialize coefficient for the first 15 Zernike terms
            Zernike_coefficient =[0,	1.00143624463888,	-51.5249359080095,	-21.9775999768718,	0.717759222249538,	15.9050351311175,	27.0075244428022,	-10.7205195257083,	-2.01254461443215,	-7.52510423390856,	10.8157068054894,	-3.02134384531165,	8.75467126995794,	-4.83968407700889,	10.0627952531622,	8.69361225696698,	-0.00569977528614984,	-0.147583443076208,	-2.30484489161088,	0.982991691361985,	-0.857744799462948,	4.92455332433157,	-10.3405736785842,	-0.702115323493294,	0.256143402391469,	2.20035061964973,	7.63689639452587,	-2.43586089370370]
                #


            # for 11X11    # initialize coefficient for the first 15 Zernike terms to get a flat image for 11X11 pupil
            Zernike_coefficient[Zernike_term] = Zernike_coefficient[Zernike_term] - 25

            for index_loop in range(0,9):
                Zernike_coefficient[Zernike_term] = Zernike_coefficient[Zernike_term] + self.Z_coef
                vars_str = str(-20+abs(index_loop * self.Z_coef))
                if -1500 <= sum(Zernike_coefficient) <= 1500:
                    Zernike_coefficientM = matlab.double(Zernike_coefficient)
                    (Z_shape, Zvol_map) = self.eng.apply_zernike(Zernike_coefficientM, nargout=2)
                    Zvol_map = [one for [one] in Zvol_map]

                    Zvol_map_Mat = matlab.double(Zvol_map)

                    self.error_code = self.eng.updateDM(self.driver_info, Zvol_map_Mat)
                    if self.error_code == 0:
                        time.sleep(0.1)
                    else:
                        return

                    for phase_shift in range(0,2):
                        # phase_shift <4 for phase shifting method
                        # phase_shift <2 for standard Zernike method
                        self.draw(self.SLM_map[phase_shift])
                        time.sleep(0.1)
                        hdu = fits.HDUList()  # initialize fits object
                        record_data = []
                        record_data = self.camera.record_to_memory(num_of_frames) / 4  # :4 to make it 14 bit

                        if record_data == []:
                            self.stop_callback()
                            return

                        hdu.append(fits.PrimaryHDU(data=record_data))
                        # other header details will come in here
                        hdu[0].header['EXPTIME'] = "%i %s" % (self.t, self.time_units.currentText())
                        if index_loop<4:
                            hdu.writeto(RM_filename +'N'+ vars_str[1:]+'_Z'+ str(Zernike_term+1) + '_' + str(phase_shift) + 'Pi.fits')
                        else:
                            hdu.writeto(RM_filename + vars_str + '_Z' + str(Zernike_term+1) + '_' + str(phase_shift) + 'Pi.fits')


        self.stop_callback()
        print('DONE!!!')
        return None

    def close_DM(self):
        """
        close DM after the response Matrix measurement
        """
        if self.driver_info != None:
            self.error_code = self.eng.closeDM(self.driver_info)
            self.eng.quit()

    def record_callback(self):
        """
        Record specified number of frames
        :return:
        """
        if not self.connected:
            return
        # check if camera had already been armed
        if self.alive and self.live_view_bool:
            try:
                self.live_view_bool = False
                self.camera.live = False  # stop loop that is producing frames
                time.sleep(0.1)
                self.record_live_thread.join()
                self.live_thread.join()
                del self.record_live_thread
                del self.live_thread
            except:
                pass
        elif not self.alive:
            self.alive = True
            self.camera.arm_camera()
            self.camera.start_recording()
            self.camera.allocate_buffer()
            self.camera._prepare_to_record_to_memory()
            self.alive = True
        else:
            pass

        hdu = fits.HDUList()  # initialize fits object

        filename = QtGui.QFileDialog.getSaveFileName(self, 'Save as..', self.save_dir)
        self.save_dir = os.path.dirname(filename)
        print(filename)

        num_of_frames = 10
        try:
            num_of_frames = int(self.FramesLab.text())
        except ValueError:
            return
        self.measurement_status.setText('Recording %d frames..'%num_of_frames)
        record_data = []
        record_data = self.camera.record_to_memory(num_of_frames)/4  # :4 to make it 14 bit

        if record_data == []:
            self.stop_callback()
            return

        self.measurement_status.setText('Exporting to FITS file')
        hdu.append(fits.PrimaryHDU(data=record_data))
        # other header details will come in here
        hdu[0].header['EXPTIME'] = "%i %s" % (self.t, self.time_units.currentText())
        hdu.writeto(filename+'.fits')
        self.measurement_status.setText('Recording finished.')
        self.stop_callback()
        return None

    def stop_callback(self):
        """
        Stops live preview
        :return:
        """
        self.alive = False

        if self.live_view_bool:
            self.live_view_bool = False
            self.camera.live = False  # stop loop that is producing frames
            time.sleep(0.1)
            self.record_live_thread.join()
            self.live_thread.join()
            del self.record_live_thread
            del self.live_thread
        # disarm camera
        self.camera.disarm_camera()
        self.display_status.setText('Idle')
        self.measurement_status.setText('')
        return

if __name__ == '__main__':

    app = QtGui.QApplication(sys.argv)
    window = QtGui.QMainWindow()
    window.setWindowTitle('PCO.PixelFly                    -ETH Zurich- ')
    try:
        icon = QtGui.QIcon('App.ico')
        window.setWindowIcon(icon)
    except:
        pass
    pco_ui = CameraWidget(parent=None)
    pco_ui.create_gui(window)
    window.show()
    sys.exit(app.exec_())







