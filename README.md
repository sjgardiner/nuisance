# NUISANCE special release for MINERvA-NOvA workshop(s) 
[Indico page, 29 Sep 2018](https://minerva-docdb.fnal.gov/cgi-bin/private/DisplayMeeting?conferenceid=6457)

----
## Description
NUISANCE release used for MINERvA-NOvA workshop on 29 Sep 2018. Used to compare the NOvA tuned GENIE 2.12.10 interaction simulation to MINERvA's tuned GENIE 2.8.6/4.

Talk is available on [FNAL DocDB](https://minerva-docdb.fnal.gov/cgi-bin/private/ShowDocument?docid=20487) and [here](http://www.hep.ph.ic.ac.uk/~cvw09/files/MINERvA_NOvA_WS_forNOVA.pdf).

Requires GENIE 2.8.4 and above and the generated GENIE files. 

----
## Instructions
Follow the usual instructions to install NUISANCE.

In a nutshell:

1. Install GENIE to some directory, export that directory as ``${GENIE}``
2. Source the ``environment_setup.sh`` file for GENIE, setting up all its dependencies
3. Get NUISANCE, get this branch, set the ``${NUISANCE}`` environment variable to the nuisance folder
4. Source ``${NUISANCE}/myenv.sh`` which sets up all your environment variables (e.g. ROOT) and generator dependencies (e.g. GENIE) — an example file is provided
5. Source ``${NUISANCE}/freshbuild.sh``, making sure it contains ``-DUSE_MINERvA_RW=1`` (for parts of the MINERvA tune), and ``-DUSE_GENIE=1`` (for GENIE event libraries and linking), with any other generators you want included
6. NUISANCE should now be built in ``${NUISANCE}/build`` with executables in ``${NUISANCE}/build/app``
7. ``${NUISANCE}/build/app/PrepareGENIE`` prepares GENIE output files for NUISANCE use. It converts the GHep format into something more manageable and produces cross-sections per target and per process. Its syntax is found by ``./PrepareGENIE --help``

----
## Running flat-tree converters, external data comparisons, and producing plots
The folder ``${NUISANCE}/MINERvA_NOvA_WS`` contains ``xml`` cardfiles for NUISANCE to produce external data comparison plots. It also contains a file ``minerva_tune.xml`` which is used to apply the MINERvA tune.

The NOvA tune is applied by the code looking for the ``nova_wgts`` branch in the GENIE file. If it can't find it it moves on, if it does it applies the total weight from that file as a reweight.

All the below assumes your raw GENIE file has been converted using ``PrepareGENIE``, outlined in step 7 above.

### To run the flat-tree converter:
Do ``./nuisflat -i "GENIE:${GENIE_PREPARED_FILE}" -o "${GENIE_PREPARED_FILE%%.root}_FlatTree_Tune.root"
`` Add ``-c minerva_tune.xml`` if you want to apply the MINERvA tune.

The flat-tree variables are intuitively named, and their meaning can be found in ``${NUISANCE}/src/MCStudies/GenericFlux_Vectors.cxx``.

To apply the tunings use the tree entries ``Weight`` for the MINERvA tune and ``CustomWeight`` for the NOvA tune. The ``CustomWeightArray[6]`` contains the different parts of the NOvA tune and are ordered as

```
fNUISANCEEvent->CustomWeight = NOVAw;
fNUISANCEEvent->CustomWeightArray[0] = MAQEw;
fNUISANCEEvent->CustomWeightArray[1] = NonResw;
fNUISANCEEvent->CustomWeightArray[2] = RPAQEw;
fNUISANCEEvent->CustomWeightArray[3] = RPARESw;
fNUISANCEEvent->CustomWeightArray[4] = MECw;
fNUISANCEEvent->CustomWeightArray[5] = NOVAw;
```
For the MINERvA tune the `CustomWeight` and `CustomWeightArray` are all 1.0.

To scale from the event rates to cross sections apply `fScaleFactor` to your distributions. All entries have their own ``fScaleFactor`` but for GENIE they are all the same, so doing ``double scale = tree->GetMinimum("fScaleFactor")`` and applying that to a given histogram provides the scaling.


### To compare against external data:
Do ``./nuiscomp -c cardfile.xml -o GENIE_FILE.root`` where the `cardfile.xml` provides NUISANCE with what experiments you want to run and where the generated events live.

The example cardfile ``${NUISANCE}/MINERvA_NOvA_WS/exp_list.xml`` provides the general structure of the `xml` cardfiles.

To find all supported experiments run ``${NUISANCE}/scripts/nuissamples``

### Producing plots:
#### Flat-tree
An example flat-tree analysis ROOT scripts is provided in `${NUISANCE}/MINERvA_NOvA_WS/FlatAnalysis.cc`. It takes two arguments:

```
root -b -q -l 'FlatAnalysis.cc("MINERvA_FlatTree.root", "NOvA_FlatTree.root")'
```

where `MINERvA_FlatTree.root` is the flat-tree produced with the MINERvA tuning. `NOvA_FlatTree.root` is the equivalent for the NOvA tuning.

#### Data comparisons
``${NUISANCE}/MINERvA_NOvA_WS/plotscript.cc`` takes four arguments: 

```
root -b -q -l 'plotscript.cc("MINERvA_untuned.root", "MINERvA_tuned.root", "NOvA_untuned.root", "NOvA_tuned.root")'
```

which are the output files from `nuiscomp` for the MINERvA untuned, tuned, NOvA untuned and tuned generator files.

These are far from perfect and exclude some known problems with some data distributions. So please contact us below if you need a hand or have requests!


----
## Known issues
The cross section scaling variable `fScaleFactor` does not play well with the NOvA target. Will be fixed.

Some of the data distributions have broken in moving beyond v2r7. E.g. MiniBooNE CCQE Q2 seems to have issues under some ROOT versions (e.g. 5.34/36). Will be fixed.

----
## Contact
Join our [Slack channel](nuisance-xsec.slack.com)

Visit our [website](nuisance.hepforge.org)

Send an email to [us](mailto:nuisance@projects.hepforge.org), [Clarence](mailto:cwret@fnal.gov), [Luke](mailto:picker24@fnal.gov), and [Kevin](mailto:kevin@rochester.edu)

