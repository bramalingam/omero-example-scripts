#-----------------------------------------------------------------------------
#  Copyright (C) 2017 University of Dundee. All rights reserved.
#
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
#  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
#------------------------------------------------------------------------------

# This Jython script uses ImageJ to Subtract Background.
# The purpose of the script is to be used in the Scripting Dialog
# of Fiji.
# Error handling is omitted to ease the reading of the script but this should be added
# if used in production to make sure the services are closed
# Information can be found at https://docs.openmicroscopy.org/omero/5.4.1/developers/Java.html

import os
from os import path

from java.lang import Long
from java.lang import String
from java.lang import System
from java.lang import Math
from java.lang import Byte
from java.util import ArrayList
from java.util import Arrays
from java.lang.reflect import Array
from jarray import zeros, array
import java

# Omero Dependencies
from omero.gateway import Gateway
from omero.gateway import LoginCredentials
from omero.gateway import SecurityContext
from omero.gateway.facility import BrowseFacility
from omero.gateway.facility import AdminFacility
from omero.gateway.facility import DataManagerFacility
from omero.gateway.facility import RawDataFacility
from omero.gateway.facility import ROIFacility
from omero.gateway.model import DatasetData
from omero.log import Logger
from omero.log import SimpleLogger

from omero.model import ExperimenterGroupI
from org.openmicroscopy.shoola.util.roi.io import ROIReader

from ome.formats.importer import ImportConfig
from ome.formats.importer import OMEROWrapper
from ome.formats.importer import ImportLibrary
from ome.formats.importer import ImportCandidates
from ome.formats.importer.cli import ErrorHandler
from ome.formats.importer.cli import LoggingImportMonitor
import loci.common
from loci.formats.in import DefaultMetadataOptions
from loci.formats.in import MetadataLevel
from loci.formats import FormatTools, ImageTools
from loci.common import DataTools
from loci.plugins.util import ImageProcessorReader

from ij import IJ, ImagePlus, ImageStack, CompositeImage
from ij.process import ByteProcessor, ShortProcessor
from ij.process import ImageProcessor


# Setup
# =====

# OMERO Server details
HOST = "outreach.openmicroscopy.org"
PORT = 4064
group_id = -1
#  parameters to edit
image_id = 4528
USERNAME = ""
PASSWORD = ""

# If you want to do analysis for someone else:
# Mention their username here (you need to have light-admin previliges)
target_user = "user-1"

# Connection method: returns a gateway object
def connect_to_omero():
    "Connect to OMERO"

    credentials = LoginCredentials()
    credentials.getServer().setHostname(HOST)
    credentials.getServer().setPort(PORT)
    credentials.getUser().setUsername(USERNAME.strip())
    credentials.getUser().setPassword(PASSWORD.strip())
    simpleLogger = SimpleLogger()
    gateway = Gateway(simpleLogger)
    
    user = gateway.connect(credentials)
    print user.getGroupId()
    return gateway

# Get a list of image_ids retured for the dataset_id
def get_image_ids(gateway, ctx, dataset_id):
    "List all image's ids contained in a Dataset"

    browse = gateway.getFacility(BrowseFacility)
    ctx = switchSecurityContext(ctx, target_user)
    
    ids = ArrayList(1)
    val = Long(dataset_id)
    ids.add(val)
    images = browse.getImagesForDatasets(ctx, ids)

    j = images.iterator()
    image_ids = []
    while j.hasNext():
        image = j.next()
        image_ids.append(String.valueOf(image.getId()))
    return image_ids


# Switch context to target user
def switchSecurityContext(ctx, target_user):
    user = gateway.getFacility(AdminFacility).lookupExperimenter(ctx, target_user)
    ctx = SecurityContext(user.getGroupId())
    ctx.setExperimenter(user)
    return ctx

#Convert omero Image object as ImageJ ImagePlus object (An alternative to OmeroReader)
def openOmeroImage(ctx, image_id):
    browse = gateway.getFacility(BrowseFacility)
    image = browse.getImage(ctx, long(image_id))
    pixels = image.getDefaultPixels()
    sizeZ = pixels.getSizeZ()
    sizeT = pixels.getSizeT()
    sizeC = pixels.getSizeC()
    sizeX = pixels.getSizeX()
    sizeY = pixels.getSizeY()
    pixtype = pixels.getPixelType()
    pixType = FormatTools.pixelTypeFromString(pixtype)
    bpp = FormatTools.getBytesPerPixel(pixType)
    isSigned = FormatTools.isSigned(pixType)
    isFloat = FormatTools.isFloatingPoint(pixType)
    isLittle = False
    interleave = False
    
    store = gateway.getPixelsStore(ctx)
    pixelsId = pixels.getId()
    store.setPixelsId(pixelsId, False)
    stack = ImageStack(sizeX, sizeY)
    for t in range(0,sizeT):
        for z in range(0,sizeZ):
            for c in range(0, sizeC):
                plane = store.getPlane(z, c, t)

                channel = ImageTools.splitChannels(plane, 0, 1, bpp, False, interleave)
                pixels = DataTools.makeDataArray(plane, bpp, isFloat, isLittle)

                q = pixels
                if (len(plane) != sizeX*sizeY):
                    tmp = q
                    q = zeros(sizeX*sizeY, 'h')
                    System.arraycopy(tmp, 0, q, 0, Math.min(len(q), len(tmp)))
                    if isSigned:
                        q = DataTools.makeSigned(q)
                    
                if q.typecode == 'b':
                    ip = ByteProcessor(sizeX, sizeY, q, None)
                elif q.typecode == 'h':
                    ip = ShortProcessor(sizeX, sizeY, q, None)
                stack.addSlice('', ip)
    # Do something
    imp = ImagePlus(image.getName(), stack)
    imp.setDimensions(sizeC, sizeZ, sizeT)
    imp.setOpenAsHyperStack(True)
    imp.show()
    return imp

# upload image to OMERO server
def upload_image(path, gateway, id):
    "Upload an image to omero"

    user = gateway.getLoggedInUser()
    sessionKey = gateway.getSessionId(user)

    config = ImportConfig()
    config.debug.set('false')
    config.hostname.set(HOST)
    config.sessionKey.set(sessionKey)
    value = "omero.model.Dataset:"
    value += str(id)
    config.target.set(value)

    loci.common.DebugTools.enableLogging("DEBUG")

    store = config.createStore()
    reader = OMEROWrapper(config)

    library = ImportLibrary(store, reader)
    error_handler = ErrorHandler(config)

    library.addObserver(LoggingImportMonitor())
    candidates = ImportCandidates(reader, path, error_handler)
    reader.setMetadataOptions(DefaultMetadataOptions(MetadataLevel.ALL))
    return library.importCandidates(config, candidates)


# Prototype analysis example
gateway = connect_to_omero()
ctx = SecurityContext(group_id)
exp = gateway.getLoggedInUser()
exp_id = exp.getId()

# get all images_ids in an omero dataset
dataset_id = 1107
ids = get_image_ids(gateway, ctx, dataset_id)

#if target_user ~= None:
# Switch context to target user and open omeroImage as ImagePlus object
ctx = switchSecurityContext(ctx, target_user)
imp = openOmeroImage(ctx, image_id)

#Some analysis which creates ROI's and Results Table
IJ.setAutoThreshold(imp, "Default dark")
IJ.run(imp, "Analyze Particles...", "clear add stack")

#Save ROI's back to OMERO
reader = ROIReader()
roi_list = reader.readImageJROIFromSources(image_id, imp)
roi_facility = gateway.getFacility(ROIFacility)
result = roi_facility.saveROIs(ctx, image_id, exp_id, roi_list)
