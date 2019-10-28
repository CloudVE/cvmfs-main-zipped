=============================
Galaxy wrapper for HiFive
=============================

This wrapper allows use of the hifive package, a set of tools for handling HiC
and 5C data, on Galaxy. This includes managing data from mapped reads, either
in bam, mat, or raw formats. All stages use hdf5 dictionaries for fast access
and minimal memory and storage usage.

This toolset includes methods for normalizing data from either HiC or 5C
experiments at the fragment-end, or fragment level resolution, respectively.
Once normalized, data can be used for plotting, binning, or other statistical
tests within the package very quickly.

Original documentation can be found `here <http://bxlab-hifive.readthedocs.org/en/latest/>`_.

Automated Installation
======================

Galaxy should be able to automatically install the dependencies, i.e. the
``h5py`` and ``hifive`` tool depedency packages.

Manual Installation
===================

For those not using Galaxy's automated installation from the Tool Shed, put
the XML and Python files in the ``tools/hifive/`` folder and add the
XML files to your ``tool_conf.xml`` as normal.  For example, use::

  <section name="HiFive" id="hifive_tools">
    <tool file="hifive.xml" />
  </section>


================
Wrapper Licence
================

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
