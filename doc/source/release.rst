.. _release:

.. currentmodule:: pandas

.. ipython:: python
   :suppress:

   import pandas as pd
   import numpy as np
   np.random.seed(123456)
   np.set_printoptions(precision=4, suppress=True)
   import matplotlib.pyplot as plt
   plt.close('all')

   pd.options.display.max_rows=15
   import pandas.util.testing as tm

*************
Release Notes
*************

This is the list of changes to pandas between each release. For full details,
see the commit logs at http://github.com/pandas-dev/pandas

**What is it**

pandas is a Python package providing fast, flexible, and expressive data
structures designed to make working with “relational” or “labeled” data both
easy and intuitive. It aims to be the fundamental high-level building block for
doing practical, real world data analysis in Python. Additionally, it has the
broader goal of becoming the most powerful and flexible open source data
analysis / manipulation tool available in any language.

**Where to get it**

* Source code: http://github.com/pandas-dev/pandas
* Binary installers on PyPI: https://pypi.org/project/pandas
* Documentation: http://pandas.pydata.org

pandas 0.23.2
-------------

**Release date**: July 5, 2018

This is a minor bug-fix release in the 0.23.x series and includes some small regression fixes
and bug fixes.

See the :ref:`full whatsnew <whatsnew_0232>` for a list of all the changes.

Thanks
~~~~~~

A total of 17 people contributed to this release.  People with a "+" by their
names contributed a patch for the first time.

* David Krych
* Jacopo Rota +
* Jeff Reback
* Jeremy Schendel
* Joris Van den Bossche
* Kalyan Gokhale
* Matthew Roeschke
* Michael Odintsov +
* Ming Li
* Pietro Battiston
* Tom Augspurger
* Uddeshya Singh
* Vu Le +
* alimcmaster1 +
* david-liu-brattle-1 +
* gfyoung
* jbrockmendel

pandas 0.23.1
-------------

**Release date**: June 12, 2018

This is a minor release from 0.23.0 and includes a number of bug fixes and
performance improvements.

See the :ref:`full whatsnew <whatsnew_0231>` for a list of all the changes.

Thanks
~~~~~~

A total of 30 people contributed to this release.  People with a "+" by their
names contributed a patch for the first time.

* Adam J. Stewart
* Adam Kim +
* Aly Sivji
* Chalmer Lowe +
* Damini Satya +
* Dr. Irv
* Gabe Fernando +
* Giftlin Rajaiah
* Jeff Reback
* Jeremy Schendel +
* Joris Van den Bossche
* Kalyan Gokhale +
* Kevin Sheppard
* Matthew Roeschke
* Max Kanter +
* Ming Li
* Pyry Kovanen +
* Stefano Cianciulli
* Tom Augspurger
* Uddeshya Singh +
* Wenhuan
* William Ayd
* chris-b1
* gfyoung
* h-vetinari
* nprad +
* ssikdar1 +
* tmnhat2001
* topper-123
* zertrin +

pandas 0.23.0
-------------

**Release date**: May 15, 2018

This is a major release from 0.22.0 and includes a number of API changes, new
features, enhancements, and performance improvements along with a large number
of bug fixes. We recommend that all users upgrade to this version.

Highlights include:

- :ref:`Round-trippable JSON format with 'table' orient <whatsnew_0230.enhancements.round-trippable_json>`.
- :ref:`Instantiation from dicts respects order for Python 3.6+ <whatsnew_0230.api_breaking.dict_insertion_order>`.
- :ref:`Dependent column arguments for assign <whatsnew_0230.enhancements.assign_dependent>`.
- :ref:`Merging / sorting on a combination of columns and index levels <whatsnew_0230.enhancements.merge_on_columns_and_levels>`.
- :ref:`Extending Pandas with custom types <whatsnew_023.enhancements.extension>`.
- :ref:`Excluding unobserved categories from groupby <whatsnew_0230.enhancements.categorical_grouping>`.
- :ref:`Changes to make output shape of DataFrame.apply consistent <whatsnew_0230.api_breaking.apply>`.

See the :ref:`full whatsnew <whatsnew_0230>` for a list of all the changes.

Thanks
~~~~~~

A total of 328 people contributed to this release.  People with a "+" by their
names contributed a patch for the first time.

* Aaron Critchley
* AbdealiJK +
* Adam Hooper +
* Albert Villanova del Moral
* Alejandro Giacometti +
* Alejandro Hohmann +
* Alex Rychyk
* Alexander Buchkovsky
* Alexander Lenail +
* Alexander Michael Schade
* Aly Sivji +
* Andreas Költringer +
* Andrew
* Andrew Bui +
* András Novoszáth +
* Andy Craze +
* Andy R. Terrel
* Anh Le +
* Anil Kumar Pallekonda +
* Antoine Pitrou +
* Antonio Linde +
* Antonio Molina +
* Antonio Quinonez +
* Armin Varshokar +
* Artem Bogachev +
* Avi Sen +
* Azeez Oluwafemi +
* Ben Auffarth +
* Bernhard Thiel +
* Bhavesh Poddar +
* BielStela +
* Blair +
* Bob Haffner
* Brett Naul +
* Brock Mendel
* Bryce Guinta +
* Carlos Eduardo Moreira dos Santos +
* Carlos García Márquez +
* Carol Willing
* Cheuk Ting Ho +
* Chitrank Dixit +
* Chris
* Chris Burr +
* Chris Catalfo +
* Chris Mazzullo
* Christian Chwala +
* Cihan Ceyhan +
* Clemens Brunner
* Colin +
* Cornelius Riemenschneider
* Crystal Gong +
* DaanVanHauwermeiren
* Dan Dixey +
* Daniel Frank +
* Daniel Garrido +
* Daniel Sakuma +
* DataOmbudsman +
* Dave Hirschfeld
* Dave Lewis +
* David Adrián Cañones Castellano +
* David Arcos +
* David C Hall +
* David Fischer
* David Hoese +
* David Lutz +
* David Polo +
* David Stansby
* Dennis Kamau +
* Dillon Niederhut
* Dimitri +
* Dr. Irv
* Dror Atariah
* Eric Chea +
* Eric Kisslinger
* Eric O. LEBIGOT (EOL) +
* FAN-GOD +
* Fabian Retkowski +
* Fer Sar +
* Gabriel de Maeztu +
* Gianpaolo Macario +
* Giftlin Rajaiah
* Gilberto Olimpio +
* Gina +
* Gjelt +
* Graham Inggs +
* Grant Roch
* Grant Smith +
* Grzegorz Konefał +
* Guilherme Beltramini
* HagaiHargil +
* Hamish Pitkeathly +
* Hammad Mashkoor +
* Hannah Ferchland +
* Hans
* Haochen Wu +
* Hissashi Rocha +
* Iain Barr +
* Ibrahim Sharaf ElDen +
* Ignasi Fosch +
* Igor Conrado Alves de Lima +
* Igor Shelvinskyi +
* Imanflow +
* Ingolf Becker
* Israel Saeta Pérez
* Iva Koevska +
* Jakub Nowacki +
* Jan F-F +
* Jan Koch +
* Jan Werkmann
* Janelle Zoutkamp +
* Jason Bandlow +
* Jaume Bonet +
* Jay Alammar +
* Jeff Reback
* JennaVergeynst
* Jimmy Woo +
* Jing Qiang Goh +
* Joachim Wagner +
* Joan Martin Miralles +
* Joel Nothman
* Joeun Park +
* John Cant +
* Johnny Metz +
* Jon Mease
* Jonas Schulze +
* Jongwony +
* Jordi Contestí +
* Joris Van den Bossche
* José F. R. Fonseca +
* Jovixe +
* Julio Martinez +
* Jörg Döpfert
* KOBAYASHI Ittoku +
* Kate Surta +
* Kenneth +
* Kevin Kuhl
* Kevin Sheppard
* Krzysztof Chomski
* Ksenia +
* Ksenia Bobrova +
* Kunal Gosar +
* Kurtis Kerstein +
* Kyle Barron +
* Laksh Arora +
* Laurens Geffert +
* Leif Walsh
* Liam Marshall +
* Liam3851 +
* Licht Takeuchi
* Liudmila +
* Ludovico Russo +
* Mabel Villalba +
* Manan Pal Singh +
* Manraj Singh
* Marc +
* Marc Garcia
* Marco Hemken +
* Maria del Mar Bibiloni +
* Mario Corchero +
* Mark Woodbridge +
* Martin Journois +
* Mason Gallo +
* Matias Heikkilä +
* Matt Braymer-Hayes
* Matt Kirk +
* Matt Maybeno +
* Matthew Kirk +
* Matthew Rocklin +
* Matthew Roeschke
* Matthias Bussonnier +
* Max Mikhaylov +
* Maxim Veksler +
* Maximilian Roos
* Maximiliano Greco +
* Michael Penkov
* Michael Röttger +
* Michael Selik +
* Michael Waskom
* Mie~~~
* Mike Kutzma +
* Ming Li +
* Mitar +
* Mitch Negus +
* Montana Low +
* Moritz Münst +
* Mortada Mehyar
* Myles Braithwaite +
* Nate Yoder
* Nicholas Ursa +
* Nick Chmura
* Nikos Karagiannakis +
* Nipun Sadvilkar +
* Nis Martensen +
* Noah +
* Noémi Éltető +
* Olivier Bilodeau +
* Ondrej Kokes +
* Onno Eberhard +
* Paul Ganssle +
* Paul Mannino +
* Paul Reidy
* Paulo Roberto de Oliveira Castro +
* Pepe Flores +
* Peter Hoffmann
* Phil Ngo +
* Pietro Battiston
* Pranav Suri +
* Priyanka Ojha +
* Pulkit Maloo +
* README Bot +
* Ray Bell +
* Riccardo Magliocchetti +
* Ridhwan Luthra +
* Robert Meyer
* Robin
* Robin Kiplang'at +
* Rohan Pandit +
* Rok Mihevc +
* Rouz Azari
* Ryszard T. Kaleta +
* Sam Cohan
* Sam Foo
* Samir Musali +
* Samuel Sinayoko +
* Sangwoong Yoon
* SarahJessica +
* Sharad Vijalapuram +
* Shubham Chaudhary +
* SiYoungOh +
* Sietse Brouwer
* Simone Basso +
* Stefania Delprete +
* Stefano Cianciulli +
* Stephen Childs +
* StephenVoland +
* Stijn Van Hoey +
* Sven
* Talitha Pumar +
* Tarbo Fukazawa +
* Ted Petrou +
* Thomas A Caswell
* Tim Hoffmann +
* Tim Swast
* Tom Augspurger
* Tommy +
* Tulio Casagrande +
* Tushar Gupta +
* Tushar Mittal +
* Upkar Lidder +
* Victor Villas +
* Vince W +
* Vinícius Figueiredo +
* Vipin Kumar +
* WBare
* Wenhuan +
* Wes Turner
* William Ayd
* Wilson Lin +
* Xbar
* Yaroslav Halchenko
* Yee Mey
* Yeongseon Choe +
* Yian +
* Yimeng Zhang
* ZhuBaohe +
* Zihao Zhao +
* adatasetaday +
* akielbowicz +
* akosel +
* alinde1 +
* amuta +
* bolkedebruin
* cbertinato
* cgohlke
* charlie0389 +
* chris-b1
* csfarkas +
* dajcs +
* deflatSOCO +
* derestle-htwg
* discort
* dmanikowski-reef +
* donK23 +
* elrubio +
* fivemok +
* fjdiod
* fjetter +
* froessler +
* gabrielclow
* gfyoung
* ghasemnaddaf
* h-vetinari +
* himanshu awasthi +
* ignamv +
* jayfoad +
* jazzmuesli +
* jbrockmendel
* jen w +
* jjames34 +
* joaoavf +
* joders +
* jschendel
* juan huguet +
* l736x +
* luzpaz +
* mdeboc +
* miguelmorin +
* miker985
* miquelcamprodon +
* orereta +
* ottiP +
* peterpanmj +
* rafarui +
* raph-m +
* readyready15728 +
* rmihael +
* samghelms +
* scriptomation +
* sfoo +
* stefansimik +
* stonebig
* tmnhat2001 +
* tomneep +
* topper-123
* tv3141 +
* verakai +
* xpvpc +
* zhanghui +

pandas 0.22.0
-------------

**Release date:** December 29, 2017

This is a major release from 0.21.1 and includes a single, API-breaking change.
We recommend that all users upgrade to this version after carefully reading the
release note.

The only changes are:

- The sum of an empty or all-*NA* ``Series`` is now ``0``
- The product of an empty or all-*NA* ``Series`` is now ``1``
- We've added a ``min_count`` parameter to ``.sum()`` and ``.prod()`` controlling
  the minimum number of valid values for the result to be valid. If fewer than
  ``min_count`` non-*NA* values are present, the result is *NA*. The default is
  ``0``. To return ``NaN``, the 0.21 behavior, use ``min_count=1``.

See the :ref:`v0.22.0 Whatsnew <whatsnew_0220>` overview for further explanation
of all the places in the library this affects.

pandas 0.21.1
-------------

**Release date:** December 12, 2017

This is a minor bug-fix release in the 0.21.x series and includes some small
regression fixes, bug fixes and performance improvements. We recommend that all
users upgrade to this version.

Highlights include:

- Temporarily restore matplotlib datetime plotting functionality. This should
  resolve issues for users who relied implicitly on pandas to plot datetimes
  with matplotlib. See :ref:`here <whatsnew_0211.converters>`.
- Improvements to the Parquet IO functions introduced in 0.21.0. See
  :ref:`here <whatsnew_0211.enhancements.parquet>`.

See the :ref:`v0.21.1 Whatsnew <whatsnew_0211>` overview for an extensive list
of all the changes for 0.21.1.

Thanks
~~~~~~

A total of 46 people contributed to this release.  People with a "+" by their
names contributed a patch for the first time.

Contributors
============

* Aaron Critchley +
* Alex Rychyk
* Alexander Buchkovsky +
* Alexander Michael Schade +
* Chris Mazzullo
* Cornelius Riemenschneider +
* Dave Hirschfeld +
* David Fischer +
* David Stansby +
* Dror Atariah +
* Eric Kisslinger +
* Hans +
* Ingolf Becker +
* Jan Werkmann +
* Jeff Reback
* Joris Van den Bossche
* Jörg Döpfert +
* Kevin Kuhl +
* Krzysztof Chomski +
* Leif Walsh
* Licht Takeuchi
* Manraj Singh +
* Matt Braymer-Hayes +
* Michael Waskom +
* Mie~~~ +
* Peter Hoffmann +
* Robert Meyer +
* Sam Cohan +
* Sietse Brouwer +
* Sven +
* Tim Swast
* Tom Augspurger
* Wes Turner
* William Ayd +
* Yee Mey +
* bolkedebruin +
* cgohlke
* derestle-htwg +
* fjdiod +
* gabrielclow +
* gfyoung
* ghasemnaddaf +
* jbrockmendel
* jschendel
* miker985 +
* topper-123

pandas 0.21.0
-------------

**Release date:** October 27, 2017

This is a major release from 0.20.3 and includes a number of API changes,
deprecations, new features, enhancements, and performance improvements along
with a large number of bug fixes. We recommend that all users upgrade to this
version.

Highlights include:

- Integration with `Apache Parquet <https://parquet.apache.org/>`__, including a new top-level :func:`read_parquet` function and :meth:`DataFrame.to_parquet` method, see :ref:`here <whatsnew_0210.enhancements.parquet>`.
- New user-facing :class:`pandas.api.types.CategoricalDtype` for specifying
  categoricals independent of the data, see :ref:`here <whatsnew_0210.enhancements.categorical_dtype>`.
- The behavior of ``sum`` and ``prod`` on all-NaN Series/DataFrames is now consistent and no longer depends on whether `bottleneck <http://berkeleyanalytics.com/bottleneck>`__ is installed, and ``sum`` and ``prod`` on empty Series now return NaN instead of 0, see :ref:`here <whatsnew_0210.api_breaking.bottleneck>`.
- Compatibility fixes for pypy, see :ref:`here <whatsnew_0210.pypy>`.
- Additions to the ``drop``, ``reindex`` and ``rename`` API to make them more consistent, see :ref:`here <whatsnew_0210.enhancements.drop_api>`.
- Addition of the new methods ``DataFrame.infer_objects`` (see :ref:`here <whatsnew_0210.enhancements.infer_objects>`) and ``GroupBy.pipe`` (see :ref:`here <whatsnew_0210.enhancements.GroupBy_pipe>`).
- Indexing with a list of labels, where one or more of the labels is missing, is deprecated and will raise a KeyError in a future version, see :ref:`here <whatsnew_0210.api_breaking.loc>`.

See the :ref:`v0.21.0 Whatsnew <whatsnew_0210>` overview for an extensive list
of all enhancements and bugs that have been fixed in 0.21.0

Thanks
~~~~~~

A total of 206 people contributed to this release.  People with a "+" by their
names contributed a patch for the first time.

Contributors
============

* 3553x +
* Aaron Barber
* Adam Gleave +
* Adam Smith +
* AdamShamlian +
* Adrian Liaw +
* Alan Velasco +
* Alan Yee +
* Alex B +
* Alex Lubbock +
* Alex Marchenko +
* Alex Rychyk +
* Amol K +
* Andreas Winkler
* Andrew +
* Andrew 亮
* André Jonasson +
* Becky Sweger
* Berkay +
* Bob Haffner +
* Bran Yang
* Brian Tu +
* Brock Mendel +
* Carol Willing +
* Carter Green +
* Chankey Pathak +
* Chris
* Chris Billington
* Chris Filo Gorgolewski +
* Chris Kerr
* Chris M +
* Chris Mazzullo +
* Christian Prinoth
* Christian Stade-Schuldt
* Christoph Moehl +
* DSM
* Daniel Chen +
* Daniel Grady
* Daniel Himmelstein
* Dave Willmer
* David Cook
* David Gwynne
* David Read +
* Dillon Niederhut +
* Douglas Rudd
* Eric Stein +
* Eric Wieser +
* Erik Fredriksen
* Florian Wilhelm +
* Floris Kint +
* Forbidden Donut
* Gabe F +
* Giftlin +
* Giftlin Rajaiah +
* Giulio Pepe +
* Guilherme Beltramini
* Guillem Borrell +
* Hanmin Qin +
* Hendrik Makait +
* Hugues Valois
* Hussain Tamboli +
* Iva Miholic +
* Jan Novotný +
* Jan Rudolph
* Jean Helie +
* Jean-Baptiste Schiratti +
* Jean-Mathieu Deschenes
* Jeff Knupp +
* Jeff Reback
* Jeff Tratner
* JennaVergeynst
* JimStearns206
* Joel Nothman
* John W. O'Brien
* Jon Crall +
* Jon Mease
* Jonathan J. Helmus +
* Joris Van den Bossche
* JosephWagner
* Juarez Bochi
* Julian Kuhlmann +
* Karel De Brabandere
* Kassandra Keeton +
* Keiron Pizzey +
* Keith Webber
* Kernc
* Kevin Sheppard
* Kirk Hansen +
* Licht Takeuchi +
* Lucas Kushner +
* Mahdi Ben Jelloul +
* Makarov Andrey +
* Malgorzata Turzanska +
* Marc Garcia +
* Margaret Sy +
* MarsGuy +
* Matt Bark +
* Matthew Roeschke
* Matti Picus
* Mehmet Ali "Mali" Akmanalp
* Michael Gasvoda +
* Michael Penkov +
* Milo +
* Morgan Stuart +
* Morgan243 +
* Nathan Ford +
* Nick Eubank
* Nick Garvey +
* Oleg Shteynbuk +
* P-Tillmann +
* Pankaj Pandey
* Patrick Luo
* Patrick O'Melveny
* Paul Reidy +
* Paula +
* Peter Quackenbush
* Peter Yanovich +
* Phillip Cloud
* Pierre Haessig
* Pietro Battiston
* Pradyumna Reddy Chinthala
* Prasanjit Prakash
* RobinFiveWords
* Ryan Hendrickson
* Sam Foo
* Sangwoong Yoon +
* Simon Gibbons +
* SimonBaron
* Steven Cutting +
* Sudeep +
* Sylvia +
* T N +
* Telt
* Thomas A Caswell
* Tim Swast +
* Tom Augspurger
* Tong SHEN
* Tuan +
* Utkarsh Upadhyay +
* Vincent La +
* Vivek +
* WANG Aiyong
* WBare
* Wes McKinney
* XF +
* Yi Liu +
* Yosuke Nakabayashi +
* aaron315 +
* abarber4gh +
* aernlund +
* agustín méndez +
* andymaheshw +
* ante328 +
* aviolov +
* bpraggastis
* cbertinato +
* cclauss +
* chernrick
* chris-b1
* dkamm +
* dwkenefick
* economy
* faic +
* fding253 +
* gfyoung
* guygoldberg +
* hhuuggoo +
* huashuai +
* ian
* iulia +
* jaredsnyder
* jbrockmendel +
* jdeschenes
* jebob +
* jschendel +
* keitakurita
* kernc +
* kiwirob +
* kjford
* linebp
* lloydkirk
* louispotok +
* majiang +
* manikbhandari +
* matthiashuschle +
* mattip
* maxwasserman +
* mjlove12 +
* nmartensen +
* pandas-docs-bot +
* parchd-1 +
* philipphanemann +
* rdk1024 +
* reidy-p +
* ri938
* ruiann +
* rvernica +
* s-weigand +
* scotthavard92 +
* skwbc +
* step4me +
* tobycheese +
* topper-123 +
* tsdlovell
* ysau +
* zzgao +


pandas 0.20.0 / 0.20.1
----------------------

**Release date:** May 5, 2017


This is a major release from 0.19.2 and includes a number of API changes, deprecations, new features,
enhancements, and performance improvements along with a large number of bug fixes. We recommend that all
users upgrade to this version.

Highlights include:

- New ``.agg()`` API for Series/DataFrame similar to the groupby-rolling-resample API's, see :ref:`here <whatsnew_0200.enhancements.agg>`
- Integration with the ``feather-format``, including a new top-level ``pd.read_feather()`` and ``DataFrame.to_feather()`` method, see :ref:`here <io.feather>`.
- The ``.ix`` indexer has been deprecated, see :ref:`here <whatsnew_0200.api_breaking.deprecate_ix>`
- ``Panel`` has been deprecated, see :ref:`here <whatsnew_0200.api_breaking.deprecate_panel>`
- Addition of an ``IntervalIndex`` and ``Interval`` scalar type, see :ref:`here <whatsnew_0200.enhancements.intervalindex>`
- Improved user API when grouping by index levels in ``.groupby()``, see :ref:`here <whatsnew_0200.enhancements.groupby_access>`
- Improved support for ``UInt64`` dtypes, see :ref:`here <whatsnew_0200.enhancements.uint64_support>`
- A new orient for JSON serialization, ``orient='table'``, that uses the Table Schema spec and that gives the possibility for a more interactive repr in the Jupyter Notebook, see :ref:`here <whatsnew_0200.enhancements.table_schema>`
- Experimental support for exporting styled DataFrames (``DataFrame.style``) to Excel, see :ref:`here <whatsnew_0200.enhancements.style_excel>`
- Window binary corr/cov operations now return a MultiIndexed ``DataFrame`` rather than a ``Panel``, as ``Panel`` is now deprecated, see :ref:`here <whatsnew_0200.api_breaking.rolling_pairwise>`
- Support for S3 handling now uses ``s3fs``, see :ref:`here <whatsnew_0200.api_breaking.s3>`
- Google BigQuery support now uses the ``pandas-gbq`` library, see :ref:`here <whatsnew_0200.api_breaking.gbq>`

See the :ref:`v0.20.1 Whatsnew <whatsnew_0200>` overview for an extensive list
of all enhancements and bugs that have been fixed in 0.20.1.


.. note::

   This is a combined release for 0.20.0 and 0.20.1.
   Version 0.20.1 contains one additional change for backwards-compatibility with downstream projects using pandas' ``utils`` routines. (:issue:`16250`)

Thanks
~~~~~~

- abaldenko
- Adam J. Stewart
- Adrian
- adrian-stepien
- Ajay Saxena
- Akash Tandon
- Albert Villanova del Moral
- Aleksey Bilogur
- alexandercbooth
- Alexis Mignon
- Amol Kahat
- Andreas Winkler
- Andrew Kittredge
- Anthonios Partheniou
- Arco Bast
- Ashish Singal
- atbd
- bastewart
- Baurzhan Muftakhidinov
- Ben Kandel
- Ben Thayer
- Ben Welsh
- Bill Chambers
- bmagnusson
- Brandon M. Burroughs
- Brian
- Brian McFee
- carlosdanielcsantos
- Carlos Souza
- chaimdemulder
- Chris
- chris-b1
- Chris Ham
- Christopher C. Aycock
- Christoph Gohlke
- Christoph Paulik
- Chris Warth
- Clemens Brunner
- DaanVanHauwermeiren
- Daniel Himmelstein
- Dave Willmer
- David Cook
- David Gwynne
- David Hoffman
- David Krych
- dickreuter
- Diego Fernandez
- Dimitris Spathis
- discort
- Dmitry L
- Dody Suria Wijaya
- Dominik Stanczak
- Dr-Irv
- Dr. Irv
- dr-leo
- D.S. McNeil
- dubourg
- dwkenefick
- Elliott Sales de Andrade
- Ennemoser Christoph
- Francesc Alted
- Fumito Hamamura
- funnycrab
- gfyoung
- Giacomo Ferroni
- goldenbull
- Graham R. Jeffries
- Greg Williams
- Guilherme Beltramini
- Guilherme Samora
- Hao Wu
- Harshit Patni
- hesham.shabana@hotmail.com
- Ilya V. Schurov
- Iván Vallés Pérez
- Jackie Leng
- Jaehoon Hwang
- James Draper
- James Goppert
- James McBride
- James Santucci
- Jan Schulz
- Jeff Carey
- Jeff Reback
- JennaVergeynst
- Jim
- Jim Crist
- Joe Jevnik
- Joel Nothman
- John
- John Tucker
- John W. O'Brien
- John Zwinck
- jojomdt
- Jonathan de Bruin
- Jonathan Whitmore
- Jon Mease
- Jon M. Mease
- Joost Kranendonk
- Joris Van den Bossche
- Joshua Bradt
- Julian Santander
- Julien Marrec
- Jun Kim
- Justin Solinsky
- Kacawi
- Kamal Kamalaldin
- Kerby Shedden
- Kernc
- Keshav Ramaswamy
- Kevin Sheppard
- Kyle Kelley
- Larry Ren
- Leon Yin
- linebp
- Line Pedersen
- Lorenzo Cestaro
- Luca Scarabello
- Lukasz
- Mahmoud Lababidi
- manu
- manuels
- Mark Mandel
- Matthew Brett
- Matthew Roeschke
- mattip
- Matti Picus
- Matt Roeschke
- maxalbert
- Maximilian Roos
- mcocdawc
- Michael Charlton
- Michael Felt
- Michael Lamparski
- Michiel Stock
- Mikolaj Chwalisz
- Min RK
- Miroslav Šedivý
- Mykola Golubyev
- Nate Yoder
- Nathalie Rud
- Nicholas Ver Halen
- Nick Chmura
- Nolan Nichols
- nuffe
- Pankaj Pandey
- paul-mannino
- Pawel Kordek
- pbreach
- Pete Huang
- Peter
- Peter Csizsek
- Petio Petrov
- Phil Ruffwind
- Pietro Battiston
- Piotr Chromiec
- Prasanjit Prakash
- Robert Bradshaw
- Rob Forgione
- Robin
- Rodolfo Fernandez
- Roger Thomas
- Rouz Azari
- Sahil Dua
- sakkemo
- Sam Foo
- Sami Salonen
- Sarah Bird
- Sarma Tangirala
- scls19fr
- Scott Sanderson
- Sebastian Bank
- Sebastian Gsänger
- Sébastien de Menten
- Shawn Heide
- Shyam Saladi
- sinhrks
- Sinhrks
- Stephen Rauch
- stijnvanhoey
- Tara Adiseshan
- themrmax
- the-nose-knows
- Thiago Serafim
- Thoralf Gutierrez
- Thrasibule
- Tobias Gustafsson
- Tom Augspurger
- tomrod
- Tong Shen
- Tong SHEN
- TrigonaMinima
- tzinckgraf
- Uwe
- wandersoncferreira
- watercrossing
- wcwagner
- Wes Turner
- Wiktor Tomczak
- WillAyd
- xgdgsc
- Yaroslav Halchenko
- Yimeng Zhang
- yui-knk


pandas 0.19.2
-------------

**Release date:** December 24, 2016

This is a minor bug-fix release in the 0.19.x series and includes some small regression fixes,
bug fixes and performance improvements.

Highlights include:

- Compatibility with Python 3.6
- Added a `Pandas Cheat Sheet <https://github.com/pandas-dev/pandas/tree/master/doc/cheatsheet/Pandas_Cheat_Sheet.pdf>`__. (:issue:`13202`).

See the :ref:`v0.19.2 Whatsnew <whatsnew_0192>` page for an overview of all
bugs that have been fixed in 0.19.2.

Thanks
~~~~~~

- Ajay Saxena
- Ben Kandel
- Chris
- Chris Ham
- Christopher C. Aycock
- Daniel Himmelstein
- Dave Willmer
- Dr-Irv
- gfyoung
- hesham shabana
- Jeff Carey
- Jeff Reback
- Joe Jevnik
- Joris Van den Bossche
- Julian Santander
- Kerby Shedden
- Keshav Ramaswamy
- Kevin Sheppard
- Luca Scarabello
- Matti Picus
- Matt Roeschke
- Maximilian Roos
- Mykola Golubyev
- Nate Yoder
- Nicholas Ver Halen
- Pawel Kordek
- Pietro Battiston
- Rodolfo Fernandez
- sinhrks
- Tara Adiseshan
- Tom Augspurger
- wandersoncferreira
- Yaroslav Halchenko


pandas 0.19.1
-------------

**Release date:** November 3, 2016

This is a minor bug-fix release from 0.19.0 and includes some small regression fixes,
bug fixes and performance improvements.

See the :ref:`v0.19.1 Whatsnew <whatsnew_0191>` page for an overview of all
bugs that have been fixed in 0.19.1.

Thanks
~~~~~~

- Adam Chainz
- Anthonios Partheniou
- Arash Rouhani
- Ben Kandel
- Brandon M. Burroughs
- Chris
- chris-b1
- Chris Warth
- David Krych
- dubourg
- gfyoung
- Iván Vallés Pérez
- Jeff Reback
- Joe Jevnik
- Jon M. Mease
- Joris Van den Bossche
- Josh Owen
- Keshav Ramaswamy
- Larry Ren
- mattrijk
- Michael Felt
- paul-mannino
- Piotr Chromiec
- Robert Bradshaw
- Sinhrks
- Thiago Serafim
- Tom Bird


pandas 0.19.0
-------------

**Release date:** October 2, 2016

This is a major release from 0.18.1 and includes number of API changes, several new features,
enhancements, and performance improvements along with a large number of bug fixes. We recommend that all
users upgrade to this version.

Highlights include:

- :func:`merge_asof` for asof-style time-series joining, see :ref:`here <whatsnew_0190.enhancements.asof_merge>`
- ``.rolling()`` is now time-series aware, see :ref:`here <whatsnew_0190.enhancements.rolling_ts>`
- :func:`read_csv` now supports parsing ``Categorical`` data, see :ref:`here <whatsnew_0190.enhancements.read_csv_categorical>`
- A function :func:`union_categorical` has been added for combining categoricals, see :ref:`here <whatsnew_0190.enhancements.union_categoricals>`
- ``PeriodIndex`` now has its own ``period`` dtype, and changed to be more consistent with other ``Index`` classes. See :ref:`here <whatsnew_0190.api.period>`
- Sparse data structures gained enhanced support of ``int`` and ``bool`` dtypes, see :ref:`here <whatsnew_0190.sparse>`
- Comparison operations with ``Series`` no longer ignores the index, see :ref:`here <whatsnew_0190.api.series_ops>` for an overview of the API changes.
- Introduction of a pandas development API for utility functions, see :ref:`here <whatsnew_0190.dev_api>`.
- Deprecation of ``Panel4D`` and ``PanelND``. We recommend to represent these types of n-dimensional data with the `xarray package <http://xarray.pydata.org/en/stable/>`__.
- Removal of the previously deprecated modules ``pandas.io.data``, ``pandas.io.wb``, ``pandas.tools.rplot``.

See the :ref:`v0.19.0 Whatsnew <whatsnew_0190>` overview for an extensive list
of all enhancements and bugs that have been fixed in 0.19.0.

Thanks
~~~~~~

- adneu
- Adrien Emery
- agraboso
- Alex Alekseyev
- Alex Vig
- Allen Riddell
- Amol
- Amol Agrawal
- Andy R. Terrel
- Anthonios Partheniou
- babakkeyvani
- Ben Kandel
- Bob Baxley
- Brett Rosen
- c123w
- Camilo Cota
- Chris
- chris-b1
- Chris Grinolds
- Christian Hudon
- Christopher C. Aycock
- Chris Warth
- cmazzullo
- conquistador1492
- cr3
- Daniel Siladji
- Douglas McNeil
- Drewrey Lupton
- dsm054
- Eduardo Blancas Reyes
- Elliot Marsden
- Evan Wright
- Felix Marczinowski
- Francis T. O'Donovan
- Gábor Lipták
- Geraint Duck
- gfyoung
- Giacomo Ferroni
- Grant Roch
- Haleemur Ali
- harshul1610
- Hassan Shamim
- iamsimha
- Iulius Curt
- Ivan Nazarov
- jackieleng
- Jeff Reback
- Jeffrey Gerard
- Jenn Olsen
- Jim Crist
- Joe Jevnik
- John Evans
- John Freeman
- John Liekezer
- Johnny Gill
- John W. O'Brien
- John Zwinck
- Jordan Erenrich
- Joris Van den Bossche
- Josh Howes
- Jozef Brandys
- Kamil Sindi
- Ka Wo Chen
- Kerby Shedden
- Kernc
- Kevin Sheppard
- Matthieu Brucher
- Maximilian Roos
- Michael Scherer
- Mike Graham
- Mortada Mehyar
- mpuels
- Muhammad Haseeb Tariq
- Nate George
- Neil Parley
- Nicolas Bonnotte
- OXPHOS
- Pan Deng / Zora
- Paul
- Pauli Virtanen
- Paul Mestemaker
- Pawel Kordek
- Pietro Battiston
- pijucha
- Piotr Jucha
- priyankjain
- Ravi Kumar Nimmi
- Robert Gieseke
- Robert Kern
- Roger Thomas
- Roy Keyes
- Russell Smith
- Sahil Dua
- Sanjiv Lobo
- Sašo Stanovnik
- Shawn Heide
- sinhrks
- Sinhrks
- Stephen Kappel
- Steve Choi
- Stewart Henderson
- Sudarshan Konge
- Thomas A Caswell
- Tom Augspurger
- Tom Bird
- Uwe Hoffmann
- wcwagner
- WillAyd
- Xiang Zhang
- Yadunandan
- Yaroslav Halchenko
- YG-Riku
- Yuichiro Kaneko
- yui-knk
- zhangjinjie
- znmean
- 颜发才（Yan Facai）

pandas 0.18.1
-------------

**Release date:** (May 3, 2016)

This is a minor release from 0.18.0 and includes a large number of bug fixes
along with several new features, enhancements, and performance improvements.

Highlights include:

- ``.groupby(...)`` has been enhanced to provide convenient syntax when working with ``.rolling(..)``, ``.expanding(..)`` and ``.resample(..)`` per group, see :ref:`here <whatsnew_0181.deferred_ops>`
- ``pd.to_datetime()`` has gained the ability to assemble dates from a ``DataFrame``, see :ref:`here <whatsnew_0181.enhancements.assembling>`
- Method chaining improvements, see :ref:`here <whatsnew_0181.enhancements.method_chain>`.
- Custom business hour offset, see :ref:`here <whatsnew_0181.enhancements.custombusinesshour>`.
- Many bug fixes in the handling of ``sparse``, see :ref:`here <whatsnew_0181.sparse>`
- Expanded the :ref:`Tutorials section <tutorial-modern>` with a feature on modern pandas, courtesy of `@TomAugsburger <https://twitter.com/TomAugspurger>`__. (:issue:`13045`).

See the :ref:`v0.18.1 Whatsnew <whatsnew_0181>` overview for an extensive list
of all enhancements and bugs that have been fixed in 0.18.1.

Thanks
~~~~~~

- Andrew Fiore-Gartland
- Bastiaan
- Benoît Vinot
- Brandon Rhodes
- DaCoEx
- Drew Fustin
- Ernesto Freitas
- Filip Ter
- Gregory Livschitz
- Gábor Lipták
- Hassan Kibirige
- Iblis Lin
- Israel Saeta Pérez
- Jason Wolosonovich
- Jeff Reback
- Joe Jevnik
- Joris Van den Bossche
- Joshua Storck
- Ka Wo Chen
- Kerby Shedden
- Kieran O'Mahony
- Leif Walsh
- Mahmoud Lababidi
- Maoyuan Liu
- Mark Roth
- Matt Wittmann
- MaxU
- Maximilian Roos
- Michael Droettboom
- Nick Eubank
- Nicolas Bonnotte
- OXPHOS
- Pauli Virtanen
- Peter Waller
- Pietro Battiston
- Prabhjot Singh
- Robin Wilson
- Roger Thomas
- Sebastian Bank
- Stephen Hoover
- Tim Hopper
- Tom Augspurger
- WANG Aiyong
- Wes Turner
- Winand
- Xbar
- Yan Facai
- adneu
- ajenkins-cargometrics
- behzad nouri
- chinskiy
- gfyoung
- jeps-journal
- jonaslb
- kotrfa
- nileracecrew
- onesandzeroes
- rs2
- sinhrks
- tsdlovell

pandas 0.18.0
-------------

**Release date:** (March 13, 2016)

This is a major release from 0.17.1 and includes a small number of API changes, several new features,
enhancements, and performance improvements along with a large number of bug fixes. We recommend that all
users upgrade to this version.

Highlights include:

- Moving and expanding window functions are now methods on Series and DataFrame,
  similar to ``.groupby``, see :ref:`here <whatsnew_0180.enhancements.moments>`.
- Adding support for a ``RangeIndex`` as a specialized form of the ``Int64Index``
  for memory savings, see :ref:`here <whatsnew_0180.enhancements.rangeindex>`.
- API breaking change to the ``.resample`` method to make it more ``.groupby``
  like, see :ref:`here <whatsnew_0180.breaking.resample>`.
- Removal of support for positional indexing with floats, which was deprecated
  since 0.14.0. This will now raise a ``TypeError``, see :ref:`here <whatsnew_0180.float_indexers>`.
- The ``.to_xarray()`` function has been added for compatibility with the
  `xarray package <http://xarray.pydata.org/en/stable/>`__, see :ref:`here <whatsnew_0180.enhancements.xarray>`.
- The ``read_sas`` function has been enhanced to read ``sas7bdat`` files, see :ref:`here <whatsnew_0180.enhancements.sas>`.
- Addition of the :ref:`.str.extractall() method <whatsnew_0180.enhancements.extract>`,
  and API changes to the :ref:`.str.extract() method <whatsnew_0180.enhancements.extract>`
  and :ref:`.str.cat() method <whatsnew_0180.enhancements.strcat>`.
- ``pd.test()`` top-level nose test runner is available (:issue:`4327`).

See the :ref:`v0.18.0 Whatsnew <whatsnew_0180>` overview for an extensive list
of all enhancements and bugs that have been fixed in 0.18.0.

Thanks
~~~~~~

- ARF
- Alex Alekseyev
- Andrew McPherson
- Andrew Rosenfeld
- Anthonios Partheniou
- Anton I. Sipos
- Ben
- Ben North
- Bran Yang
- Chris
- Chris Carroux
- Christopher C. Aycock
- Christopher Scanlin
- Cody
- Da Wang
- Daniel Grady
- Dorozhko Anton
- Dr-Irv
- Erik M. Bray
- Evan Wright
- Francis T. O'Donovan
- Frank Cleary
- Gianluca Rossi
- Graham Jeffries
- Guillaume Horel
- Henry Hammond
- Isaac Schwabacher
- Jean-Mathieu Deschenes
- Jeff Reback
- Joe Jevnik
- John Freeman
- John Fremlin
- Jonas Hoersch
- Joris Van den Bossche
- Joris Vankerschaver
- Justin Lecher
- Justin Lin
- Ka Wo Chen
- Keming Zhang
- Kerby Shedden
- Kyle
- Marco Farrugia
- MasonGallo
- MattRijk
- Matthew Lurie
- Maximilian Roos
- Mayank Asthana
- Mortada Mehyar
- Moussa Taifi
- Navreet Gill
- Nicolas Bonnotte
- Paul Reiners
- Philip Gura
- Pietro Battiston
- RahulHP
- Randy Carnevale
- Rinoc Johnson
- Rishipuri
- Sangmin Park
- Scott E Lasley
- Sereger13
- Shannon Wang
- Skipper Seabold
- Thierry Moisan
- Thomas A Caswell
- Toby Dylan Hocking
- Tom Augspurger
- Travis
- Trent Hauck
- Tux1
- Varun
- Wes McKinney
- Will Thompson
- Yoav Ram
- Yoong Kang Lim
- Yoshiki Vázquez Baeza
- Young Joong Kim
- Younggun Kim
- Yuval Langer
- alex argunov
- behzad nouri
- boombard
- brian-pantano
- chromy
- daniel
- dgram0
- gfyoung
- hack-c
- hcontrast
- jfoo
- kaustuv deolal
- llllllllll
- ranarag
- rockg
- scls19fr
- seales
- sinhrks
- srib
- surveymedia.ca
- tworec

pandas 0.17.1
-------------

**Release date:** (November 21, 2015)

This is a minor release from 0.17.0 and includes a large number of bug fixes
along with several new features, enhancements, and performance improvements.

Highlights include:

- Support for Conditional HTML Formatting, see :ref:`here <whatsnew_0171.style>`
- Releasing the GIL on the csv reader & other ops, see :ref:`here <whatsnew_0171.performance>`
- Regression in ``DataFrame.drop_duplicates`` from 0.16.2, causing incorrect results on integer values (:issue:`11376`)

See the :ref:`v0.17.1 Whatsnew <whatsnew_0171>` overview for an extensive list
of all enhancements and bugs that have been fixed in 0.17.1.

Thanks
~~~~~~

- Aleksandr Drozd
- Alex Chase
- Anthonios Partheniou
- BrenBarn
- Brian J. McGuirk
- Chris
- Christian Berendt
- Christian Perez
- Cody Piersall
- Data & Code Expert Experimenting with Code on Data
- DrIrv
- Evan Wright
- Guillaume Gay
- Hamed Saljooghinejad
- Iblis Lin
- Jake VanderPlas
- Jan Schulz
- Jean-Mathieu Deschenes
- Jeff Reback
- Jimmy Callin
- Joris Van den Bossche
- K.-Michael Aye
- Ka Wo Chen
- Loïc Séguin-C
- Luo Yicheng
- Magnus Jöud
- Manuel Leonhardt
- Matthew Gilbert
- Maximilian Roos
- Michael
- Nicholas Stahl
- Nicolas Bonnotte
- Pastafarianist
- Petra Chong
- Phil Schaf
- Philipp A
- Rob deCarvalho
- Roman Khomenko
- Rémy Léone
- Sebastian Bank
- Thierry Moisan
- Tom Augspurger
- Tux1
- Varun
- Wieland Hoffmann
- Winterflower
- Yoav Ram
- Younggun Kim
- Zeke
- ajcr
- azuranski
- behzad nouri
- cel4
- emilydolson
- hironow
- lexual
- llllllllll
- rockg
- silentquasar
- sinhrks
- taeold

pandas 0.17.0
-------------

**Release date:** (October 9, 2015)

This is a major release from 0.16.2 and includes a small number of API changes, several new features,
enhancements, and performance improvements along with a large number of bug fixes. We recommend that all
users upgrade to this version.

Highlights include:

- Release the Global Interpreter Lock (GIL) on some cython operations, see :ref:`here <whatsnew_0170.gil>`
- Plotting methods are now available as attributes of the ``.plot`` accessor, see :ref:`here <whatsnew_0170.plot>`
- The sorting API has been revamped to remove some long-time inconsistencies, see :ref:`here <whatsnew_0170.api_breaking.sorting>`
- Support for a ``datetime64[ns]`` with timezones as a first-class dtype, see :ref:`here <whatsnew_0170.tz>`
- The default for ``to_datetime`` will now be to ``raise`` when presented with unparseable formats,
  previously this would return the original input. Also, date parse
  functions now return consistent results. See :ref:`here <whatsnew_0170.api_breaking.to_datetime>`
- The default for ``dropna`` in ``HDFStore`` has changed to ``False``, to store by default all rows even
  if they are all ``NaN``, see :ref:`here <whatsnew_0170.api_breaking.hdf_dropna>`
- Datetime accessor (``dt``) now supports ``Series.dt.strftime`` to generate formatted strings for datetime-likes, and ``Series.dt.total_seconds`` to generate each duration of the timedelta in seconds. See :ref:`here <whatsnew_0170.strftime>`
- ``Period`` and ``PeriodIndex`` can handle multiplied freq like ``3D``, which corresponding to 3 days span. See :ref:`here <whatsnew_0170.periodfreq>`
- Development installed versions of pandas will now have ``PEP440`` compliant version strings (:issue:`9518`)
- Development support for benchmarking with the `Air Speed Velocity library <https://github.com/spacetelescope/asv/>`_ (:issue:`8316`)
- Support for reading SAS xport files, see :ref:`here <whatsnew_0170.enhancements.sas_xport>`
- Documentation comparing SAS to *pandas*, see :ref:`here <compare_with_sas>`
- Removal of the automatic TimeSeries broadcasting, deprecated since 0.8.0, see :ref:`here <whatsnew_0170.prior_deprecations>`
- Display format with plain text can optionally align with Unicode East Asian Width, see :ref:`here <whatsnew_0170.east_asian_width>`
- Compatibility with Python 3.5 (:issue:`11097`)
- Compatibility with matplotlib 1.5.0 (:issue:`11111`)

See the :ref:`v0.17.0 Whatsnew <whatsnew_0170>` overview for an extensive list
of all enhancements and bugs that have been fixed in 0.17.0.

Thanks
~~~~~~

- Alex Rothberg
- Andrea Bedini
- Andrew Rosenfeld
- Andy Li
- Anthonios Partheniou
- Artemy Kolchinsky
- Bernard Willers
- Charlie Clark
- Chris
- Chris Whelan
- Christoph Gohlke
- Christopher Whelan
- Clark Fitzgerald
- Clearfield Christopher
- Dan Ringwalt
- Daniel Ni
- Data & Code Expert Experimenting with Code on Data
- David Cottrell
- David John Gagne
- David Kelly
- ETF
- Eduardo Schettino
- Egor
- Egor Panfilov
- Evan Wright
- Frank Pinter
- Gabriel Araujo
- Garrett-R
- Gianluca Rossi
- Guillaume Gay
- Guillaume Poulin
- Harsh Nisar
- Ian Henriksen
- Ian Hoegen
- Jaidev Deshpande
- Jan Rudolph
- Jan Schulz
- Jason Swails
- Jeff Reback
- Jonas Buyl
- Joris Van den Bossche
- Joris Vankerschaver
- Josh Levy-Kramer
- Julien Danjou
- Ka Wo Chen
- Karrie Kehoe
- Kelsey Jordahl
- Kerby Shedden
- Kevin Sheppard
- Lars Buitinck
- Leif Johnson
- Luis Ortiz
- Mac
- Matt Gambogi
- Matt Savoie
- Matthew Gilbert
- Maximilian Roos
- Michelangelo D'Agostino
- Mortada Mehyar
- Nick Eubank
- Nipun Batra
- Ondřej Čertík
- Phillip Cloud
- Pratap Vardhan
- Rafal Skolasinski
- Richard Lewis
- Rinoc Johnson
- Rob Levy
- Robert Gieseke
- Safia Abdalla
- Samuel Denny
- Saumitra Shahapure
- Sebastian Pölsterl
- Sebastian Rubbert
- Sheppard, Kevin
- Sinhrks
- Siu Kwan Lam
- Skipper Seabold
- Spencer Carrucciu
- Stephan Hoyer
- Stephen Hoover
- Stephen Pascoe
- Terry Santegoeds
- Thomas Grainger
- Tjerk Santegoeds
- Tom Augspurger
- Vincent Davis
- Winterflower
- Yaroslav Halchenko
- Yuan Tang (Terry)
- agijsberts
- ajcr
- behzad nouri
- cel4
- cyrusmaher
- davidovitch
- ganego
- jreback
- juricast
- larvian
- maximilianr
- msund
- rekcahpassyla
- robertzk
- scls19fr
- seth-p
- sinhrks
- springcoil
- terrytangyuan
- tzinckgraf

pandas 0.16.2
-------------

**Release date:** (June 12, 2015)

This is a minor release from 0.16.1 and includes a large number of bug fixes
along with several new features, enhancements, and performance improvements.

Highlights include:

- A new ``pipe`` method, see :ref:`here <whatsnew_0162.enhancements.pipe>`
- Documentation on how to use `numba <http://numba.pydata.org>`_ with *pandas*, see :ref:`here <enhancingperf.numba>`

See the :ref:`v0.16.2 Whatsnew <whatsnew_0162>` overview for an extensive list
of all enhancements and bugs that have been fixed in 0.16.2.

Thanks
~~~~~~

- Andrew Rosenfeld
- Artemy Kolchinsky
- Bernard Willers
- Christer van der Meeren
- Christian Hudon
- Constantine Glen Evans
- Daniel Julius Lasiman
- Evan Wright
- Francesco Brundu
- Gaëtan de Menten
- Jake VanderPlas
- James Hiebert
- Jeff Reback
- Joris Van den Bossche
- Justin Lecher
- Ka Wo Chen
- Kevin Sheppard
- Mortada Mehyar
- Morton Fox
- Robin Wilson
- Thomas Grainger
- Tom Ajamian
- Tom Augspurger
- Yoshiki Vázquez Baeza
- Younggun Kim
- austinc
- behzad nouri
- jreback
- lexual
- rekcahpassyla
- scls19fr
- sinhrks

pandas 0.16.1
-------------

**Release date:** (May 11, 2015)

This is a minor release from 0.16.0 and includes a large number of bug fixes
along with several new features, enhancements, and performance improvements.
A small number of API changes were necessary to fix existing bugs.

See the :ref:`v0.16.1 Whatsnew <whatsnew_0161>` overview for an extensive list
of all API changes, enhancements and bugs that have been fixed in 0.16.1.

Thanks
~~~~~~

- Alfonso MHC
- Andy Hayden
- Artemy Kolchinsky
- Chris Gilmer
- Chris Grinolds
- Dan Birken
- David BROCHART
- David Hirschfeld
- David Stephens
- Dr. Leo
- Evan Wright
- Frans van Dunné
- Hatem Nassrat
- Henning Sperr
- Hugo Herter
- Jan Schulz
- Jeff Blackburne
- Jeff Reback
- Jim Crist
- Jonas Abernot
- Joris Van den Bossche
- Kerby Shedden
- Leo Razoumov
- Manuel Riel
- Mortada Mehyar
- Nick Burns
- Nick Eubank
- Olivier Grisel
- Phillip Cloud
- Pietro Battiston
- Roy Hyunjin Han
- Sam Zhang
- Scott Sanderson
- Stephan Hoyer
- Tiago Antao
- Tom Ajamian
- Tom Augspurger
- Tomaz Berisa
- Vikram Shirgur
- Vladimir Filimonov
- William Hogman
- Yasin A
- Younggun Kim
- behzad nouri
- dsm054
- floydsoft
- flying-sheep
- gfr
- jnmclarty
- jreback
- ksanghai
- lucas
- mschmohl
- ptype
- rockg
- scls19fr
- sinhrks


pandas 0.16.0
-------------

**Release date:** (March 22, 2015)

This is a major release from 0.15.2 and includes a number of API changes, several new features, enhancements, and
performance improvements along with a large number of bug fixes.

Highlights include:

- ``DataFrame.assign`` method, see :ref:`here <whatsnew_0160.enhancements.assign>`
- ``Series.to_coo/from_coo`` methods to interact with ``scipy.sparse``, see :ref:`here <whatsnew_0160.enhancements.sparse>`
- Backwards incompatible change to ``Timedelta`` to conform the ``.seconds`` attribute with ``datetime.timedelta``, see :ref:`here <whatsnew_0160.api_breaking.timedelta>`
- Changes to the ``.loc`` slicing API to conform with the behavior of ``.ix`` see :ref:`here <whatsnew_0160.api_breaking.indexing>`
- Changes to the default for ordering in the ``Categorical`` constructor, see :ref:`here <whatsnew_0160.api_breaking.categorical>`
- The ``pandas.tools.rplot``, ``pandas.sandbox.qtpandas`` and ``pandas.rpy``
  modules are deprecated. We refer users to external packages like
  `seaborn <http://stanford.edu/~mwaskom/software/seaborn/>`_,
  `pandas-qt <https://github.com/datalyze-solutions/pandas-qt>`_ and
  `rpy2 <http://rpy2.bitbucket.org/>`_ for similar or equivalent
  functionality, see :ref:`here <whatsnew_0160.deprecations>`

See the :ref:`v0.16.0 Whatsnew <whatsnew_0160>` overview or the issue tracker on GitHub for an extensive list
of all API changes, enhancements and bugs that have been fixed in 0.16.0.

Thanks
~~~~~~

- Aaron Toth
- Alan Du
- Alessandro Amici
- Artemy Kolchinsky
- Ashwini Chaudhary
- Ben Schiller
- Bill Letson
- Brandon Bradley
- Chau Hoang
- Chris Reynolds
- Chris Whelan
- Christer van der Meeren
- David Cottrell
- David Stephens
- Ehsan Azarnasab
- Garrett-R
- Guillaume Gay
- Jake Torcasso
- Jason Sexauer
- Jeff Reback
- John McNamara
- Joris Van den Bossche
- Joschka zur Jacobsmühlen
- Juarez Bochi
- Junya Hayashi
- K.-Michael Aye
- Kerby Shedden
- Kevin Sheppard
- Kieran O'Mahony
- Kodi Arfer
- Matti Airas
- Min RK
- Mortada Mehyar
- Robert
- Scott E Lasley
- Scott Lasley
- Sergio Pascual
- Skipper Seabold
- Stephan Hoyer
- Thomas Grainger
- Tom Augspurger
- TomAugspurger
- Vladimir Filimonov
- Vyomkesh Tripathi
- Will Holmgren
- Yulong Yang
- behzad nouri
- bertrandhaut
- bjonen
- cel4
- clham
- hsperr
- ischwabacher
- jnmclarty
- josham
- jreback
- omtinez
- roch
- sinhrks
- unutbu

pandas 0.15.2
-------------

**Release date:** (December 12, 2014)

This is a minor release from 0.15.1 and includes a large number of bug fixes
along with several new features, enhancements, and performance improvements.
A small number of API changes were necessary to fix existing bugs.

See the :ref:`v0.15.2 Whatsnew <whatsnew_0152>` overview for an extensive list
of all API changes, enhancements and bugs that have been fixed in 0.15.2.

Thanks
~~~~~~

- Aaron Staple
- Angelos Evripiotis
- Artemy Kolchinsky
- Benoit Pointet
- Brian Jacobowski
- Charalampos Papaloizou
- Chris Warth
- David Stephens
- Fabio Zanini
- Francesc Via
- Henry Kleynhans
- Jake VanderPlas
- Jan Schulz
- Jeff Reback
- Jeff Tratner
- Joris Van den Bossche
- Kevin Sheppard
- Matt Suggit
- Matthew Brett
- Phillip Cloud
- Rupert Thompson
- Scott E Lasley
- Stephan Hoyer
- Stephen Simmons
- Sylvain Corlay
- Thomas Grainger
- Tiago Antao
- Trent Hauck
- Victor Chaves
- Victor Salgado
- Vikram Bhandoh
- WANG Aiyong
- Will Holmgren
- behzad nouri
- broessli
- charalampos papaloizou
- immerrr
- jnmclarty
- jreback
- mgilbert
- onesandzeroes
- peadarcoyle
- rockg
- seth-p
- sinhrks
- unutbu
- wavedatalab
- Åsmund Hjulstad

pandas 0.15.1
-------------

**Release date:** (November 9, 2014)

This is a minor release from 0.15.0 and includes a small number of API changes, several new features, enhancements, and
performance improvements along with a large number of bug fixes.

See the :ref:`v0.15.1 Whatsnew <whatsnew_0151>` overview for an extensive list
of all API changes, enhancements and bugs that have been fixed in 0.15.1.

Thanks
~~~~~~

- Aaron Staple
- Andrew Rosenfeld
- Anton I. Sipos
- Artemy Kolchinsky
- Bill Letson
- Dave Hughes
- David Stephens
- Guillaume Horel
- Jeff Reback
- Joris Van den Bossche
- Kevin Sheppard
- Nick Stahl
- Sanghee Kim
- Stephan Hoyer
- TomAugspurger
- WANG Aiyong
- behzad nouri
- immerrr
- jnmclarty
- jreback
- pallav-fdsi
- unutbu

pandas 0.15.0
-------------

**Release date:** (October 18, 2014)

This is a major release from 0.14.1 and includes a number of API changes, several new features, enhancements, and
performance improvements along with a large number of bug fixes.

Highlights include:

- Drop support for NumPy < 1.7.0 (:issue:`7711`)
- The ``Categorical`` type was integrated as a first-class pandas type, see :ref:`here <whatsnew_0150.cat>`
- New scalar type ``Timedelta``, and a new index type ``TimedeltaIndex``, see :ref:`here <whatsnew_0150.timedeltaindex>`
- New DataFrame default display for ``df.info()`` to include memory usage, see :ref:`Memory Usage <whatsnew_0150.memory>`
- New datetimelike properties accessor ``.dt`` for Series, see :ref:`Datetimelike Properties <whatsnew_0150.dt>`
- Split indexing documentation into :ref:`Indexing and Selecting Data <indexing>` and :ref:`MultiIndex / Advanced Indexing <advanced>`
- Split out string methods documentation into :ref:`Working with Text Data <text>`
- ``read_csv`` will now by default ignore blank lines when parsing, see :ref:`here <whatsnew_0150.blanklines>`
- API change in using Indexes in set operations, see :ref:`here <whatsnew_0150.index_set_ops>`
- Internal refactoring of the ``Index`` class to no longer sub-class ``ndarray``, see :ref:`Internal Refactoring <whatsnew_0150.refactoring>`
- dropping support for ``PyTables`` less than version 3.0.0, and ``numexpr`` less than version 2.1 (:issue:`7990`)

See the :ref:`v0.15.0 Whatsnew <whatsnew_0150>` overview or the issue tracker on GitHub for an extensive list
of all API changes, enhancements and bugs that have been fixed in 0.15.0.

Thanks
~~~~~~

- Aaron Schumacher
- Adam Greenhall
- Andy Hayden
- Anthony O'Brien
- Artemy Kolchinsky
- behzad nouri
- Benedikt Sauer
- benjamin
- Benjamin Thyreau
- Ben Schiller
- bjonen
- BorisVerk
- Chris Reynolds
- Chris Stoafer
- Dav Clark
- dlovell
- DSM
- dsm054
- FragLegs
- German Gomez-Herrero
- Hsiaoming Yang
- Huan Li
- hunterowens
- Hyungtae Kim
- immerrr
- Isaac Slavitt
- ischwabacher
- Jacob Schaer
- Jacob Wasserman
- Jan Schulz
- Jeff Tratner
- Jesse Farnham
- jmorris0x0
- jnmclarty
- Joe Bradish
- Joerg Rittinger
- John W. O'Brien
- Joris Van den Bossche
- jreback
- Kevin Sheppard
- klonuo
- Kyle Meyer
- lexual
- Max Chang
- mcjcode
- Michael Mueller
- Michael W Schatzow
- Mike Kelly
- Mortada Mehyar
- mtrbean
- Nathan Sanders
- Nathan Typanski
- onesandzeroes
- Paul Masurel
- Phillip Cloud
- Pietro Battiston
- RenzoBertocchi
- rockg
- Ross Petchler
- seth-p
- Shahul Hameed
- Shashank Agarwal
- sinhrks
- someben
- stahlous
- stas-sl
- Stephan Hoyer
- thatneat
- tom-alcorn
- TomAugspurger
- Tom Augspurger
- Tony Lorenzo
- unknown
- unutbu
- Wes Turner
- Wilfred Hughes
- Yevgeniy Grechka
- Yoshiki VÃ¡zquez Baeza
- zachcp

pandas 0.14.1
-------------

**Release date:** (July 11, 2014)

This is a minor release from 0.14.0 and includes a small number of API changes, several new features, enhancements, and
performance improvements along with a large number of bug fixes.

Highlights include:

- New methods :meth:`~pandas.DataFrame.select_dtypes` to select columns
  based on the dtype and :meth:`~pandas.Series.sem` to calculate the
  standard error of the mean.
- Support for dateutil timezones (see :ref:`docs <timeseries.timezone>`).
- Support for ignoring full line comments in the :func:`~pandas.read_csv`
  text parser.
- New documentation section on :ref:`Options and Settings <options>`.
- Lots of bug fixes.

See the :ref:`v0.14.1 Whatsnew <whatsnew_0141>` overview or the issue tracker on GitHub for an extensive list
of all API changes, enhancements and bugs that have been fixed in 0.14.1.

Thanks
~~~~~~

- Andrew Rosenfeld
- Andy Hayden
- Benjamin Adams
- Benjamin M. Gross
- Brian Quistorff
- Brian Wignall
- bwignall
- clham
- Daniel Waeber
- David Bew
- David Stephens
- DSM
- dsm054
- helger
- immerrr
- Jacob Schaer
- jaimefrio
- Jan Schulz
- John David Reaver
- John W. O'Brien
- Joris Van den Bossche
- jreback
- Julien Danjou
- Kevin Sheppard
- K.-Michael Aye
- Kyle Meyer
- lexual
- Matthew Brett
- Matt Wittmann
- Michael Mueller
- Mortada Mehyar
- onesandzeroes
- Phillip Cloud
- Rob Levy
- rockg
- sanguineturtle
- Schaer, Jacob C
- seth-p
- sinhrks
- Stephan Hoyer
- Thomas Kluyver
- Todd Jennings
- TomAugspurger
- unknown
- yelite

pandas 0.14.0
-------------

**Release date:** (May 31, 2014)

This is a major release from 0.13.1 and includes a number of API changes, several new features, enhancements, and
performance improvements along with a large number of bug fixes.

Highlights include:

- Officially support Python 3.4
- SQL interfaces updated to use ``sqlalchemy``, see :ref:`here<whatsnew_0140.sql>`.
- Display interface changes, see :ref:`here<whatsnew_0140.display>`
- MultiIndexing using Slicers, see :ref:`here<whatsnew_0140.slicers>`.
- Ability to join a singly-indexed DataFrame with a MultiIndexed DataFrame, see :ref:`here <merging.join_on_mi>`
- More consistency in groupby results and more flexible groupby specifications, see :ref:`here<whatsnew_0140.groupby>`
- Holiday calendars are now supported in ``CustomBusinessDay``, see :ref:`here <timeseries.holiday>`
- Several improvements in plotting functions, including: hexbin, area and pie plots, see :ref:`here<whatsnew_0140.plotting>`.
- Performance doc section on I/O operations, see :ref:`here <io.perf>`

See the :ref:`v0.14.0 Whatsnew <whatsnew_0140>` overview or the issue tracker on GitHub for an extensive list
of all API changes, enhancements and bugs that have been fixed in 0.14.0.

Thanks
~~~~~~

- Acanthostega
- Adam Marcus
- agijsberts
- akittredge
- Alex Gaudio
- Alex Rothberg
- AllenDowney
- Andrew Rosenfeld
- Andy Hayden
- ankostis
- anomrake
- Antoine Mazières
- anton-d
- bashtage
- Benedikt Sauer
- benjamin
- Brad Buran
- bwignall
- cgohlke
- chebee7i
- Christopher Whelan
- Clark Fitzgerald
- clham
- Dale Jung
- Dan Allan
- Dan Birken
- danielballan
- Daniel Waeber
- David Jung
- David Stephens
- Douglas McNeil
- DSM
- Garrett Drapala
- Gouthaman Balaraman
- Guillaume Poulin
- hshimizu77
- hugo
- immerrr
- ischwabacher
- Jacob Howard
- Jacob Schaer
- jaimefrio
- Jason Sexauer
- Jeff Reback
- Jeffrey Starr
- Jeff Tratner
- John David Reaver
- John McNamara
- John W. O'Brien
- Jonathan Chambers
- Joris Van den Bossche
- jreback
- jsexauer
- Julia Evans
- Júlio
- Katie Atkinson
- kdiether
- Kelsey Jordahl
- Kevin Sheppard
- K.-Michael Aye
- Matthias Kuhn
- Matt Wittmann
- Max Grender-Jones
- Michael E. Gruen
- michaelws
- mikebailey
- Mike Kelly
- Nipun Batra
- Noah Spies
- ojdo
- onesandzeroes
- Patrick O'Keeffe
- phaebz
- Phillip Cloud
- Pietro Battiston
- PKEuS
- Randy Carnevale
- ribonoous
- Robert Gibboni
- rockg
- sinhrks
- Skipper Seabold
- SplashDance
- Stephan Hoyer
- Tim Cera
- Tobias Brandt
- Todd Jennings
- TomAugspurger
- Tom Augspurger
- unutbu
- westurner
- Yaroslav Halchenko
- y-p
- zach powers

pandas 0.13.1
-------------

**Release date:** (February 3, 2014)

New Features
~~~~~~~~~~~~

- Added ``date_format`` and ``datetime_format`` attribute to ``ExcelWriter``.
  (:issue:`4133`)

API Changes
~~~~~~~~~~~

- ``Series.sort`` will raise a ``ValueError`` (rather than a ``TypeError``) on sorting an
  object that is a view of another (:issue:`5856`, :issue:`5853`)
- Raise/Warn ``SettingWithCopyError`` (according to the option ``chained_assignment`` in more cases,
  when detecting chained assignment, related (:issue:`5938`, :issue:`6025`)
- DataFrame.head(0) returns self instead of empty frame (:issue:`5846`)
- ``autocorrelation_plot`` now accepts ``**kwargs``. (:issue:`5623`)
- ``convert_objects`` now accepts a ``convert_timedeltas='coerce'`` argument to allow forced dtype conversion of
  timedeltas (:issue:`5458`,:issue:`5689`)
- Add ``-NaN`` and ``-nan`` to the default set of NA values
  (:issue:`5952`).  See :ref:`NA Values <io.na_values>`.
- ``NDFrame`` now has an ``equals`` method. (:issue:`5283`)
- ``DataFrame.apply`` will use the ``reduce`` argument to determine whether a
  ``Series`` or a ``DataFrame`` should be returned when the ``DataFrame`` is
  empty (:issue:`6007`).

Experimental Features
~~~~~~~~~~~~~~~~~~~~~

Improvements to existing features
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- perf improvements in Series datetime/timedelta binary operations (:issue:`5801`)
- `option_context` context manager now available as top-level API (:issue:`5752`)
- df.info() view now display dtype info per column (:issue:`5682`)
- df.info() now honors option max_info_rows, disable null counts for large frames (:issue:`5974`)
- perf improvements in DataFrame ``count/dropna`` for ``axis=1``
- Series.str.contains now has a `regex=False` keyword which can be faster for plain (non-regex) string patterns. (:issue:`5879`)
- support ``dtypes`` property on ``Series/Panel/Panel4D``
- extend ``Panel.apply`` to allow arbitrary functions (rather than only ufuncs) (:issue:`1148`)
  allow multiple axes to be used to operate on slabs of a ``Panel``
- The ``ArrayFormatter`` for ``datetime`` and ``timedelta64`` now intelligently
  limit precision based on the values in the array (:issue:`3401`)
- ``pd.show_versions()`` is now available for convenience when reporting issues.
- perf improvements to Series.str.extract (:issue:`5944`)
- perf improvements in ``dtypes/ftypes`` methods (:issue:`5968`)
- perf improvements in indexing with object dtypes (:issue:`5968`)
- improved dtype inference for ``timedelta`` like passed to constructors (:issue:`5458`, :issue:`5689`)
- escape special characters when writing to latex (:issue: `5374`)
- perf improvements in ``DataFrame.apply`` (:issue:`6013`)
- ``pd.read_csv`` and ``pd.to_datetime`` learned a new ``infer_datetime_format`` keyword which greatly
  improves parsing perf in many cases. Thanks to @lexual for suggesting and @danbirken
  for rapidly implementing. (:issue:`5490`,:issue:`6021`)
- add ability to recognize '%p' format code (am/pm) to date parsers when the specific format
  is supplied (:issue:`5361`)
- Fix performance regression in JSON IO (:issue:`5765`)
- performance regression in Index construction from Series (:issue:`6150`)

.. _release.bug_fixes-0.13.1:

Bug Fixes
~~~~~~~~~

- Bug in ``io.wb.get_countries`` not including all countries (:issue:`6008`)
- Bug in Series replace with timestamp dict (:issue:`5797`)
- read_csv/read_table now respects the `prefix` kwarg (:issue:`5732`).
- Bug in selection with missing values via ``.ix`` from a duplicate indexed DataFrame failing (:issue:`5835`)
- Fix issue of boolean comparison on empty DataFrames (:issue:`5808`)
- Bug in isnull handling ``NaT`` in an object array (:issue:`5443`)
- Bug in ``to_datetime`` when passed a ``np.nan`` or integer datelike and a format string (:issue:`5863`)
- Bug in groupby dtype conversion with datetimelike (:issue:`5869`)
- Regression in handling of empty Series as indexers to Series  (:issue:`5877`)
- Bug in internal caching, related to (:issue:`5727`)
- Testing bug in reading JSON/msgpack from a non-filepath on windows under py3 (:issue:`5874`)
- Bug when assigning to .ix[tuple(...)] (:issue:`5896`)
- Bug in fully reindexing a Panel (:issue:`5905`)
- Bug in idxmin/max with object dtypes (:issue:`5914`)
- Bug in ``BusinessDay`` when adding n days to a date not on offset when n>5 and n%5==0 (:issue:`5890`)
- Bug in assigning to chained series with a series via ix (:issue:`5928`)
- Bug in creating an empty DataFrame, copying, then assigning (:issue:`5932`)
- Bug in DataFrame.tail with empty frame (:issue:`5846`)
- Bug in propagating metadata on ``resample`` (:issue:`5862`)
- Fixed string-representation of ``NaT`` to be "NaT" (:issue:`5708`)
- Fixed string-representation for Timestamp to show nanoseconds if present (:issue:`5912`)
- ``pd.match`` not returning passed sentinel
- ``Panel.to_frame()`` no longer fails when ``major_axis`` is a
  ``MultiIndex`` (:issue:`5402`).
- Bug in ``pd.read_msgpack`` with inferring a ``DateTimeIndex`` frequency
  incorrectly (:issue:`5947`)
- Fixed ``to_datetime`` for array with both Tz-aware datetimes and ``NaT``'s  (:issue:`5961`)
- Bug in rolling skew/kurtosis when passed a Series with bad data (:issue:`5749`)
- Bug in scipy ``interpolate`` methods with a datetime index (:issue:`5975`)
- Bug in NaT comparison if a mixed datetime/np.datetime64 with NaT were passed (:issue:`5968`)
- Fixed bug with ``pd.concat`` losing dtype information if all inputs are empty (:issue:`5742`)
- Recent changes in IPython cause warnings to be emitted when using previous versions
  of pandas in QTConsole, now fixed. If you're using an older version and
  need to suppress the warnings, see (:issue:`5922`).
- Bug in merging ``timedelta`` dtypes (:issue:`5695`)
- Bug in plotting.scatter_matrix function. Wrong alignment among diagonal
  and off-diagonal plots, see (:issue:`5497`).
- Regression in Series with a MultiIndex via ix (:issue:`6018`)
- Bug in Series.xs with a MultiIndex (:issue:`6018`)
- Bug in Series construction of mixed type with datelike and an integer (which should result in
  object type and not automatic conversion) (:issue:`6028`)
- Possible segfault when chained indexing with an object array under NumPy 1.7.1 (:issue:`6026`, :issue:`6056`)
- Bug in setting using fancy indexing a single element with a non-scalar (e.g. a list),
  (:issue:`6043`)
- ``to_sql`` did not respect ``if_exists`` (:issue:`4110` :issue:`4304`)
- Regression in ``.get(None)`` indexing from 0.12 (:issue:`5652`)
- Subtle ``iloc`` indexing bug, surfaced in (:issue:`6059`)
- Bug with insert of strings into DatetimeIndex (:issue:`5818`)
- Fixed unicode bug in to_html/HTML repr (:issue:`6098`)
- Fixed missing arg validation in get_options_data (:issue:`6105`)
- Bug in assignment with duplicate columns in a frame where the locations
  are a slice (e.g. next to each other) (:issue:`6120`)
- Bug in propagating _ref_locs during construction of a DataFrame with dups
  index/columns (:issue:`6121`)
- Bug in ``DataFrame.apply`` when using mixed datelike reductions (:issue:`6125`)
- Bug in ``DataFrame.append`` when appending a row with different columns (:issue:`6129`)
- Bug in DataFrame construction with recarray and non-ns datetime dtype (:issue:`6140`)
- Bug in ``.loc`` setitem indexing with a dataframe on rhs, multiple item setting, and
  a datetimelike (:issue:`6152`)
- Fixed a bug in ``query``/``eval`` during lexicographic string comparisons (:issue:`6155`).
- Fixed a bug in ``query`` where the index of a single-element ``Series`` was
  being thrown away (:issue:`6148`).
- Bug in ``HDFStore`` on appending a dataframe with MultiIndexed columns to
  an existing table (:issue:`6167`)
- Consistency with dtypes in setting an empty DataFrame (:issue:`6171`)
- Bug in selecting on a MultiIndex ``HDFStore`` even in the presence of under
  specified column spec (:issue:`6169`)
- Bug in ``nanops.var`` with ``ddof=1`` and 1 elements would sometimes return ``inf``
  rather than ``nan`` on some platforms (:issue:`6136`)
- Bug in Series and DataFrame bar plots ignoring the ``use_index`` keyword (:issue:`6209`)
- Bug in groupby with mixed str/int under python3 fixed; ``argsort`` was failing (:issue:`6212`)

pandas 0.13.0
-------------

**Release date:** January 3, 2014

New Features
~~~~~~~~~~~~

- ``plot(kind='kde')`` now accepts the optional parameters ``bw_method`` and
  ``ind``, passed to scipy.stats.gaussian_kde() (for scipy >= 0.11.0) to set
  the bandwidth, and to gkde.evaluate() to specify the indices at which it
  is evaluated, respectively. See scipy docs. (:issue:`4298`)
- Added ``isin`` method to DataFrame (:issue:`4211`)
- ``df.to_clipboard()`` learned a new ``excel`` keyword that let's you
  paste df data directly into excel (enabled by default). (:issue:`5070`).
- Clipboard functionality now works with PySide (:issue:`4282`)
- New ``extract`` string method returns regex matches more conveniently
  (:issue:`4685`)
- Auto-detect field widths in read_fwf when unspecified (:issue:`4488`)
- ``to_csv()`` now outputs datetime objects according to a specified format
  string via the ``date_format`` keyword (:issue:`4313`)
- Added ``LastWeekOfMonth`` DateOffset (:issue:`4637`)
- Added ``cumcount`` groupby method (:issue:`4646`)
- Added ``FY5253``, and ``FY5253Quarter`` DateOffsets (:issue:`4511`)
- Added ``mode()`` method to ``Series`` and ``DataFrame`` to get the
  statistical mode(s) of a column/series. (:issue:`5367`)

Experimental Features
~~~~~~~~~~~~~~~~~~~~~

- The new :func:`~pandas.eval` function implements expression evaluation
  using ``numexpr`` behind the scenes. This results in large speedups for
  complicated expressions involving large DataFrames/Series.
- :class:`~pandas.DataFrame` has a new :meth:`~pandas.DataFrame.eval` that
  evaluates an expression in the context of the ``DataFrame``; allows
  inline expression assignment
- A :meth:`~pandas.DataFrame.query` method has been added that allows
  you to select elements of a ``DataFrame`` using a natural query syntax
  nearly identical to Python syntax.
- ``pd.eval`` and friends now evaluate operations involving ``datetime64``
  objects in Python space because ``numexpr`` cannot handle ``NaT`` values
  (:issue:`4897`).
- Add msgpack support via ``pd.read_msgpack()`` and ``pd.to_msgpack()`` /
  ``df.to_msgpack()`` for serialization of arbitrary pandas (and python
  objects) in a lightweight portable binary format (:issue:`686`, :issue:`5506`)
- Added PySide support for the qtpandas DataFrameModel and DataFrameWidget.
- Added :mod:`pandas.io.gbq` for reading from (and writing to) Google
  BigQuery into a DataFrame. (:issue:`4140`)

Improvements to existing features
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- ``read_html`` now raises a ``URLError`` instead of catching and raising a
  ``ValueError`` (:issue:`4303`, :issue:`4305`)
- ``read_excel`` now supports an integer in its ``sheetname`` argument giving
  the index of the sheet to read in (:issue:`4301`).
- ``get_dummies`` works with NaN (:issue:`4446`)
- Added a test for ``read_clipboard()`` and ``to_clipboard()``
  (:issue:`4282`)
- Added bins argument to ``value_counts`` (:issue:`3945`), also sort and
  ascending, now available in Series method as well as top-level function.
- Text parser now treats anything that reads like inf ("inf", "Inf", "-Inf",
  "iNf", etc.) to infinity. (:issue:`4220`, :issue:`4219`), affecting
  ``read_table``, ``read_csv``, etc.
- Added a more informative error message when plot arguments contain
  overlapping color and style arguments (:issue:`4402`)
- Significant table writing performance improvements in ``HDFStore``
- JSON date serialization now performed in low-level C code.
- JSON support for encoding datetime.time
- Expanded JSON docs, more info about orient options and the use of the numpy
  param when decoding.
- Add ``drop_level`` argument to xs (:issue:`4180`)
- Can now resample a DataFrame with ohlc (:issue:`2320`)
- ``Index.copy()`` and ``MultiIndex.copy()`` now accept keyword arguments to
  change attributes (i.e., ``names``, ``levels``, ``labels``)
  (:issue:`4039`)
- Add ``rename`` and ``set_names`` methods to ``Index`` as well as
  ``set_names``, ``set_levels``, ``set_labels`` to ``MultiIndex``.
  (:issue:`4039`) with improved validation for all (:issue:`4039`,
  :issue:`4794`)
- A Series of dtype ``timedelta64[ns]`` can now be divided/multiplied
  by an integer series (:issue:`4521`)
- A Series of dtype ``timedelta64[ns]`` can now be divided by another
  ``timedelta64[ns]`` object to yield a ``float64`` dtyped Series. This
  is frequency conversion; astyping is also supported.
- Timedelta64 support ``fillna/ffill/bfill`` with an integer interpreted as
  seconds, or a ``timedelta`` (:issue:`3371`)
- Box numeric ops on ``timedelta`` Series (:issue:`4984`)
- Datetime64 support ``ffill/bfill``
- Performance improvements with ``__getitem__`` on ``DataFrames`` with
  when the key is a column
- Support for using a ``DatetimeIndex/PeriodsIndex`` directly in a datelike
  calculation e.g. s-s.index (:issue:`4629`)
- Better/cleaned up exceptions in core/common, io/excel and core/format
  (:issue:`4721`, :issue:`3954`), as well as cleaned up test cases in
  tests/test_frame, tests/test_multilevel (:issue:`4732`).
- Performance improvement of timeseries plotting with PeriodIndex and added
  test to vbench (:issue:`4705` and :issue:`4722`)
- Add ``axis`` and ``level`` keywords to ``where``, so that the ``other``
  argument can now be an alignable pandas object.
- ``to_datetime`` with a format of '%Y%m%d' now parses much faster
- It's now easier to hook new Excel writers into pandas (just subclass
  ``ExcelWriter`` and register your engine). You can specify an ``engine`` in
  ``to_excel`` or in ``ExcelWriter``.  You can also specify which writers you
  want to use by default with config options ``io.excel.xlsx.writer`` and
  ``io.excel.xls.writer``. (:issue:`4745`, :issue:`4750`)
- ``Panel.to_excel()`` now accepts keyword arguments that will be passed to
  its ``DataFrame``'s ``to_excel()`` methods. (:issue:`4750`)
- Added XlsxWriter as an optional ``ExcelWriter``  engine. This is about 5x
  faster than the default openpyxl xlsx writer and is equivalent in speed
  to the xlwt xls writer module. (:issue:`4542`)
- allow DataFrame constructor to accept more list-like objects, e.g. list of
  ``collections.Sequence`` and ``array.Array`` objects (:issue:`3783`,
  :issue:`4297`, :issue:`4851`), thanks @lgautier
- DataFrame constructor now accepts a NumPy masked record array
  (:issue:`3478`), thanks @jnothman
- ``__getitem__`` with ``tuple`` key (e.g., ``[:, 2]``) on ``Series``
  without ``MultiIndex`` raises ``ValueError`` (:issue:`4759`, :issue:`4837`)
- ``read_json`` now raises a (more informative) ``ValueError`` when the dict
  contains a bad key and ``orient='split'`` (:issue:`4730`, :issue:`4838`)
- ``read_stata`` now accepts Stata 13 format (:issue:`4291`)
- ``ExcelWriter`` and ``ExcelFile`` can be used as context managers.
  (:issue:`3441`, :issue:`4933`)
- ``pandas`` is now tested with two different versions of ``statsmodels``
  (0.4.3 and 0.5.0) (:issue:`4981`).
- Better string representations of ``MultiIndex`` (including ability to
  roundtrip via ``repr``). (:issue:`3347`, :issue:`4935`)
- Both ExcelFile and read_excel to accept an xlrd.Book for the io
  (formerly path_or_buf) argument; this requires engine to be set.
  (:issue:`4961`).
- ``concat`` now gives a more informative error message when passed objects
  that cannot be concatenated (:issue:`4608`).
- Add ``halflife`` option to exponentially weighted moving functions (PR
  :issue:`4998`)
- ``to_dict`` now takes ``records`` as a possible out type.  Returns an array
  of column-keyed dictionaries. (:issue:`4936`)
- ``tz_localize`` can infer a fall daylight savings transition based on the
  structure of unlocalized data (:issue:`4230`)
- DatetimeIndex is now in the API documentation
- Improve support for converting R datasets to pandas objects (more
  informative index for timeseries and numeric, support for factors, dist,
  and high-dimensional arrays).
- :func:`~pandas.read_html` now supports the ``parse_dates``,
  ``tupleize_cols`` and ``thousands`` parameters (:issue:`4770`).
- :meth:`~pandas.io.json.json_normalize` is a new method to allow you to
  create a flat table from semi-structured JSON data. :ref:`See the
  docs<io.json_normalize>` (:issue:`1067`)
- ``DataFrame.from_records()`` will now accept generators (:issue:`4910`)
- ``DataFrame.interpolate()`` and ``Series.interpolate()`` have been expanded
  to include interpolation methods from scipy. (:issue:`4434`, :issue:`1892`)
- ``Series`` now supports a ``to_frame`` method to convert it to a
  single-column DataFrame (:issue:`5164`)
- DatetimeIndex (and date_range) can now be constructed in a left- or
  right-open fashion using the ``closed`` parameter (:issue:`4579`)
- Python csv parser now supports usecols (:issue:`4335`)
- Added support for Google Analytics v3 API segment IDs that also supports v2
  IDs. (:issue:`5271`)
- ``NDFrame.drop()`` now accepts names as well as integers for the axis
  argument. (:issue:`5354`)
- Added short docstrings to a few methods that were missing them + fixed the
  docstrings for Panel flex methods. (:issue:`5336`)
- ``NDFrame.drop()``, ``NDFrame.dropna()``, and ``.drop_duplicates()`` all
  accept ``inplace`` as a keyword argument; however, this only means that the
  wrapper is updated inplace, a copy is still made internally.
  (:issue:`1960`, :issue:`5247`, :issue:`5628`, and related :issue:`2325` [still not
  closed])
- Fixed bug in `tools.plotting.andrews_curvres` so that lines are drawn grouped
  by color as expected.
- ``read_excel()`` now tries to convert integral floats (like ``1.0``) to int
  by default. (:issue:`5394`)
- Excel writers now have a default option ``merge_cells`` in ``to_excel()``
  to merge cells in MultiIndex and Hierarchical Rows. Note: using this
  option it is no longer possible to round trip Excel files with merged
  MultiIndex and Hierarchical Rows. Set the ``merge_cells`` to ``False`` to
  restore the previous behaviour.  (:issue:`5254`)
- The FRED DataReader now accepts multiple series (:issue:`3413`)
- StataWriter adjusts variable names to Stata's limitations (:issue:`5709`)

API Changes
~~~~~~~~~~~

- ``DataFrame.reindex()`` and forward/backward filling now raises ValueError
  if either index is not monotonic (:issue:`4483`, :issue:`4484`).
- ``pandas`` now is Python 2/3 compatible without the need for 2to3 thanks to
  @jtratner. As a result, pandas now uses iterators more extensively. This
  also led to the introduction of substantive parts of the Benjamin
  Peterson's ``six`` library into compat. (:issue:`4384`, :issue:`4375`,
  :issue:`4372`)
- ``pandas.util.compat`` and ``pandas.util.py3compat`` have been merged into
  ``pandas.compat``. ``pandas.compat`` now includes many functions allowing
  2/3 compatibility. It contains both list and iterator versions of range,
  filter, map and zip, plus other necessary elements for Python 3
  compatibility. ``lmap``, ``lzip``, ``lrange`` and ``lfilter`` all produce
  lists instead of iterators, for compatibility with ``numpy``, subscripting
  and ``pandas`` constructors.(:issue:`4384`, :issue:`4375`, :issue:`4372`)
- deprecated ``iterkv``, which will be removed in a future release (was just
  an alias of iteritems used to get around ``2to3``'s changes).
  (:issue:`4384`, :issue:`4375`, :issue:`4372`)
- ``Series.get`` with negative indexers now returns the same as ``[]``
  (:issue:`4390`)
- allow ``ix/loc`` for Series/DataFrame/Panel to set on any axis even when
  the single-key is not currently contained in the index for that axis
  (:issue:`2578`, :issue:`5226`, :issue:`5632`, :issue:`5720`,
  :issue:`5744`, :issue:`5756`)
- Default export for ``to_clipboard`` is now csv with a sep of `\t` for
  compat (:issue:`3368`)
- ``at`` now will enlarge the object inplace (and return the same)
  (:issue:`2578`)
- ``DataFrame.plot`` will scatter plot x versus y by passing
  ``kind='scatter'`` (:issue:`2215`)

- ``HDFStore``

  - ``append_to_multiple`` automatically synchronizes writing rows to multiple
    tables and adds a ``dropna`` kwarg (:issue:`4698`)
  - handle a passed ``Series`` in table format (:issue:`4330`)
  - added an ``is_open`` property to indicate if the underlying file handle
    is_open; a closed store will now report 'CLOSED' when viewing the store
    (rather than raising an error) (:issue:`4409`)
  - a close of a ``HDFStore`` now will close that instance of the
    ``HDFStore`` but will only close the actual file if the ref count (by
    ``PyTables``) w.r.t. all of the open handles are 0. Essentially you have
    a local instance of ``HDFStore`` referenced by a variable. Once you close
    it, it will report closed. Other references (to the same file) will
    continue to operate until they themselves are closed. Performing an
    action on a closed file will raise ``ClosedFileError``
  - removed the ``_quiet`` attribute, replace by a ``DuplicateWarning`` if
    retrieving duplicate rows from a table (:issue:`4367`)
  - removed the ``warn`` argument from ``open``. Instead a
    ``PossibleDataLossError`` exception will be raised if you try to use
    ``mode='w'`` with an OPEN file handle (:issue:`4367`)
  - allow a passed locations array or mask as a ``where`` condition
    (:issue:`4467`)
  - add the keyword ``dropna=True`` to ``append`` to change whether ALL nan
    rows are not written to the store (default is ``True``, ALL nan rows are
    NOT written), also settable via the option ``io.hdf.dropna_table``
    (:issue:`4625`)
  - the ``format`` keyword now replaces the ``table`` keyword; allowed values
    are ``fixed(f)|table(t)`` the ``Storer`` format has been renamed to
    ``Fixed``
  - a column MultiIndex will be recreated properly (:issue:`4710`); raise on
    trying to use a MultiIndex with data_columns on the same axis
  - ``select_as_coordinates`` will now return an ``Int64Index`` of the
    resultant selection set
  - support ``timedelta64[ns]`` as a serialization type (:issue:`3577`)
  - store `datetime.date` objects as ordinals rather then time-tuples to avoid
    timezone issues (:issue:`2852`), thanks @tavistmorph and @numpand
  - ``numexpr`` 2.2.2 fixes incompatibility in PyTables 2.4 (:issue:`4908`)
  - ``flush`` now accepts an ``fsync`` parameter, which defaults to ``False``
    (:issue:`5364`)
  - ``unicode`` indices not supported on ``table`` formats (:issue:`5386`)
  - pass through store creation arguments; can be used to support in-memory stores
- ``JSON``

  - added ``date_unit`` parameter to specify resolution of timestamps.
    Options are seconds, milliseconds, microseconds and nanoseconds.
    (:issue:`4362`, :issue:`4498`).
  - added ``default_handler`` parameter to allow a callable to be passed
    which will be responsible for handling otherwise unserialiable objects.
    (:issue:`5138`)

- ``Index`` and ``MultiIndex`` changes (:issue:`4039`):

  - Setting ``levels`` and ``labels`` directly on ``MultiIndex`` is now
    deprecated. Instead, you can use the ``set_levels()`` and
    ``set_labels()`` methods.
  - ``levels``, ``labels`` and ``names`` properties no longer return lists,
    but instead return containers that do not allow setting of items
    ('mostly immutable')
  - ``levels``, ``labels`` and ``names`` are validated upon setting and are
    either copied or shallow-copied.
  - inplace setting of ``levels`` or ``labels`` now correctly invalidates the
    cached properties. (:issue:`5238`).
  - ``__deepcopy__`` now returns a shallow copy (currently: a view) of the
    data - allowing metadata changes.
  - ``MultiIndex.astype()`` now only allows ``np.object_``-like dtypes and
    now returns a ``MultiIndex`` rather than an ``Index``. (:issue:`4039`)
  - Added ``is_`` method to ``Index`` that allows fast equality comparison of
    views (similar to ``np.may_share_memory`` but no false positives, and
    changes on ``levels`` and ``labels`` setting on ``MultiIndex``).
    (:issue:`4859` , :issue:`4909`)
  - Aliased ``__iadd__`` to ``__add__``. (:issue:`4996`)
  - Added ``is_`` method to ``Index`` that allows fast equality comparison of
    views (similar to ``np.may_share_memory`` but no false positives, and
    changes on ``levels`` and ``labels`` setting on ``MultiIndex``).
    (:issue:`4859`, :issue:`4909`)

- Infer and downcast dtype if ``downcast='infer'`` is passed to
  ``fillna/ffill/bfill`` (:issue:`4604`)
- ``__nonzero__`` for all NDFrame objects, will now raise a ``ValueError``,
  this reverts back to (:issue:`1073`, :issue:`4633`) behavior. Add
  ``.bool()`` method to ``NDFrame`` objects to facilitate evaluating of
  single-element boolean Series
- ``DataFrame.update()`` no longer raises a ``DataConflictError``, it now
  will raise a ``ValueError`` instead (if necessary) (:issue:`4732`)
- ``Series.isin()`` and ``DataFrame.isin()``  now raise a ``TypeError`` when
  passed a string (:issue:`4763`). Pass a ``list`` of one element (containing
  the string) instead.
- Remove undocumented/unused ``kind`` keyword argument from ``read_excel``,
  and ``ExcelFile``. (:issue:`4713`, :issue:`4712`)
- The ``method`` argument of ``NDFrame.replace()`` is valid again, so that a
  a list can be passed to ``to_replace`` (:issue:`4743`).
- provide automatic dtype conversions on _reduce operations (:issue:`3371`)
- exclude non-numerics if mixed types with datelike in _reduce operations
  (:issue:`3371`)
- default for ``tupleize_cols`` is now ``False`` for both ``to_csv`` and
  ``read_csv``. Fair warning in 0.12 (:issue:`3604`)
- moved timedeltas support to pandas.tseries.timedeltas.py; add timedeltas
  string parsing, add top-level ``to_timedelta`` function
- ``NDFrame`` now is compatible with Python's toplevel ``abs()`` function
  (:issue:`4821`).
- raise a ``TypeError`` on invalid comparison ops on Series/DataFrame (e.g.
  integer/datetime) (:issue:`4968`)
- Added a new index type, ``Float64Index``. This will be automatically
  created when passing floating values in index creation.  This enables a
  pure label-based slicing paradigm that makes ``[],ix,loc`` for scalar
  indexing and slicing work exactly the same.  Indexing on other index types
  are preserved (and positional fall back for ``[],ix``), with the exception,
  that floating point slicing on indexes on non ``Float64Index`` will raise a
  ``TypeError``, e.g. ``Series(range(5))[3.5:4.5]`` (:issue:`263`,:issue:`5375`)
- Make Categorical repr nicer (:issue:`4368`)
- Remove deprecated ``Factor`` (:issue:`3650`)
- Remove deprecated ``set_printoptions/reset_printoptions`` (:issue:``3046``)
- Remove deprecated ``_verbose_info`` (:issue:`3215`)
- Begin removing methods that don't make sense on ``GroupBy`` objects
  (:issue:`4887`).
- Remove deprecated ``read_clipboard/to_clipboard/ExcelFile/ExcelWriter``
  from ``pandas.io.parsers`` (:issue:`3717`)
- All non-Index NDFrames (``Series``, ``DataFrame``, ``Panel``, ``Panel4D``,
  ``SparsePanel``, etc.), now support the entire set of arithmetic operators
  and arithmetic flex methods (add, sub, mul, etc.). ``SparsePanel`` does not
  support ``pow`` or ``mod`` with non-scalars. (:issue:`3765`)
- Arithmetic func factories are now passed real names (suitable for using
  with super) (:issue:`5240`)
- Provide NumPy compatibility with 1.7 for a calling convention like
  ``np.prod(pandas_object)`` as NumPy call with additional keyword args
  (:issue:`4435`)
- Provide __dir__ method (and local context) for tab completion / remove
  ipython completers code (:issue:`4501`)
- Support non-unique axes in a Panel via indexing operations (:issue:`4960`)
- ``.truncate`` will raise a ``ValueError`` if invalid before and afters
  dates are given (:issue:`5242`)
- ``Timestamp`` now supports ``now/today/utcnow`` class methods
  (:issue:`5339`)
- default for `display.max_seq_len` is now 100 rather then `None`. This activates
  truncated display ("...") of long sequences in various places. (:issue:`3391`)
- **All** division with ``NDFrame`` - likes is now true division, regardless
  of the future import. You can use ``//`` and ``floordiv`` to do integer
  division.

.. code-block:: ipython

   In [3]: arr = np.array([1, 2, 3, 4])

   In [4]: arr2 = np.array([5, 3, 2, 1])

   In [5]: arr / arr2
   Out[5]: array([0, 0, 1, 4])

   In [6]: pd.Series(arr) / pd.Series(arr2) # no future import required
   Out[6]:
   0    0.200000
   1    0.666667
   2    1.500000
   3    4.000000
   dtype: float64

- raise/warn ``SettingWithCopyError/Warning`` exception/warning when setting of a
  copy through chained assignment is detected, settable via option ``mode.chained_assignment``
- test the list of ``NA`` values in the csv parser. add ``N/A``, ``#NA`` as independent default
  na values (:issue:`5521`)
- The refactoring involving ``Series`` deriving from ``NDFrame`` breaks ``rpy2<=2.3.8``. an Issue
  has been opened against rpy2 and a workaround is detailed in :issue:`5698`. Thanks @JanSchulz.
- ``Series.argmin`` and ``Series.argmax`` are now aliased to ``Series.idxmin`` and ``Series.idxmax``.
  These return the *index* of the min or max element respectively. Prior to 0.13.0 these would return
  the position of the min / max element (:issue:`6214`)

Internal Refactoring
~~~~~~~~~~~~~~~~~~~~

In 0.13.0 there is a major refactor primarily to subclass ``Series`` from
``NDFrame``, which is the base class currently for ``DataFrame`` and ``Panel``,
to unify methods and behaviors. Series formerly subclassed directly from
``ndarray``. (:issue:`4080`, :issue:`3862`, :issue:`816`)
See :ref:`Internal Refactoring<whatsnew_0130.refactoring>`

- Refactor of series.py/frame.py/panel.py to move common code to generic.py

 - added ``_setup_axes`` to created generic NDFrame structures
 - moved methods

   - ``from_axes``, ``_wrap_array``, ``axes``, ``ix``, ``loc``, ``iloc``,
     ``shape``, ``empty``, ``swapaxes``, ``transpose``, ``pop``
   - ``__iter__``, ``keys``, ``__contains__``, ``__len__``, ``__neg__``,
     ``__invert__``
   - ``convert_objects``, ``as_blocks``, ``as_matrix``, ``values``
   - ``__getstate__``, ``__setstate__`` (compat remains in frame/panel)
   - ``__getattr__``, ``__setattr__``
   - ``_indexed_same``, ``reindex_like``, ``align``, ``where``, ``mask``
   - ``fillna``, ``replace`` (``Series`` replace is now consistent with
     ``DataFrame``)
   - ``filter`` (also added axis argument to selectively filter on a different
     axis)
   - ``reindex``, ``reindex_axis``, ``take``
   - ``truncate`` (moved to become part of ``NDFrame``)
   - ``isnull/notnull`` now available on ``NDFrame`` objects

- These are API changes which make ``Panel`` more consistent with ``DataFrame``

 - ``swapaxes`` on a ``Panel`` with the same axes specified now return a copy
 - support attribute access for setting
 - ``filter`` supports same API as original ``DataFrame`` filter
 - ``fillna`` refactored to ``core/generic.py``, while > 3ndim is
   ``NotImplemented``

- Series now inherits from ``NDFrame`` rather than directly from ``ndarray``.
  There are several minor changes that affect the API.

 - NumPy functions that do not support the array interface will now return
   ``ndarrays`` rather than series, e.g. ``np.diff``, ``np.ones_like``,
   ``np.where``
 - ``Series(0.5)`` would previously return the scalar ``0.5``, this is no
   longer supported
 - ``TimeSeries`` is now an alias for ``Series``. the property
   ``is_time_series`` can be used to distinguish (if desired)

- Refactor of Sparse objects to use BlockManager

 - Created a new block type in internals, ``SparseBlock``, which can hold
   multi-dtypes and is non-consolidatable. ``SparseSeries`` and
   ``SparseDataFrame`` now inherit more methods from there hierarchy
   (Series/DataFrame), and no longer inherit from ``SparseArray`` (which
   instead is the object of the ``SparseBlock``)
 - Sparse suite now supports integration with non-sparse data. Non-float
   sparse data is supportable (partially implemented)
 - Operations on sparse structures within DataFrames should preserve
   sparseness, merging type operations will convert to dense (and back to
   sparse), so might be somewhat inefficient
 - enable setitem on ``SparseSeries`` for boolean/integer/slices
 - ``SparsePanels`` implementation is unchanged (e.g. not using BlockManager,
   needs work)

- added ``ftypes`` method to Series/DataFame, similar to ``dtypes``, but
  indicates if the underlying is sparse/dense (as well as the dtype)
- All ``NDFrame`` objects now have a ``_prop_attributes``, which can be used
  to indicate various values to propagate to a new object from an existing
  (e.g. name in ``Series`` will follow more automatically now)
- Internal type checking is now done via a suite of generated classes,
  allowing ``isinstance(value, klass)`` without having to directly import the
  klass, courtesy of @jtratner
- Bug in Series update where the parent frame is not updating its cache based
  on changes (:issue:`4080`, :issue:`5216`) or types (:issue:`3217`), fillna
  (:issue:`3386`)
- Indexing with dtype conversions fixed (:issue:`4463`, :issue:`4204`)
- Refactor ``Series.reindex`` to core/generic.py (:issue:`4604`,
  :issue:`4618`), allow ``method=`` in reindexing on a Series to work
- ``Series.copy`` no longer accepts the ``order`` parameter and is now
  consistent with ``NDFrame`` copy
- Refactor ``rename`` methods to core/generic.py; fixes ``Series.rename`` for
  (:issue:`4605`), and adds ``rename`` with the same signature for ``Panel``
- Series (for index) / Panel (for items) now as attribute access to its
  elements  (:issue:`1903`)
- Refactor ``clip`` methods to core/generic.py (:issue:`4798`)
- Refactor of ``_get_numeric_data/_get_bool_data`` to core/generic.py,
  allowing Series/Panel functionality
- Refactor of Series arithmetic with time-like objects
  (datetime/timedelta/time etc.) into a separate, cleaned up wrapper class.
  (:issue:`4613`)
- Complex compat for ``Series`` with ``ndarray``. (:issue:`4819`)
- Removed unnecessary ``rwproperty`` from code base in favor of builtin
  property. (:issue:`4843`)
- Refactor object level numeric methods (mean/sum/min/max...) from object
  level modules to ``core/generic.py`` (:issue:`4435`).
- Refactor cum objects to core/generic.py (:issue:`4435`), note that these
  have a more numpy-like function signature.
- :func:`~pandas.read_html` now uses ``TextParser`` to parse HTML data from
  bs4/lxml (:issue:`4770`).
- Removed the ``keep_internal`` keyword parameter in
  ``pandas/core/groupby.py`` because it wasn't being used (:issue:`5102`).
- Base ``DateOffsets`` are no longer all instantiated on importing pandas,
  instead they are generated and cached on the fly. The internal
  representation and handling of DateOffsets has also been clarified.
  (:issue:`5189`, related :issue:`5004`)
- ``MultiIndex`` constructor now validates that passed levels and labels are
  compatible. (:issue:`5213`, :issue:`5214`)
- Unity ``dropna`` for Series/DataFrame signature (:issue:`5250`),
  tests from :issue:`5234`, courtesy of @rockg
- Rewrite assert_almost_equal() in cython for performance (:issue:`4398`)
- Added an internal ``_update_inplace`` method to facilitate updating
  ``NDFrame`` wrappers on inplace ops (only is for convenience of caller,
  doesn't actually prevent copies). (:issue:`5247`)

.. _release.bug_fixes-0.13.0:


Bug Fixes
~~~~~~~~~

- ``HDFStore``

  - raising an invalid ``TypeError`` rather than ``ValueError`` when
    appending with a different block ordering (:issue:`4096`)
  - ``read_hdf`` was not respecting as passed ``mode`` (:issue:`4504`)
  - appending a 0-len table will work correctly (:issue:`4273`)
  - ``to_hdf`` was raising when passing both arguments ``append`` and
    ``table`` (:issue:`4584`)
  - reading from a store with duplicate columns across dtypes would raise
    (:issue:`4767`)
  - Fixed a bug where ``ValueError`` wasn't correctly raised when column
    names weren't strings (:issue:`4956`)
  - A zero length series written in Fixed format not deserializing properly.
    (:issue:`4708`)
  - Fixed decoding perf issue on pyt3 (:issue:`5441`)
  - Validate levels in a MultiIndex before storing (:issue:`5527`)
  - Correctly handle ``data_columns`` with a Panel (:issue:`5717`)
- Fixed bug in tslib.tz_convert(vals, tz1, tz2): it could raise IndexError
  exception while trying to access trans[pos + 1] (:issue:`4496`)
- The ``by`` argument now works correctly with the ``layout`` argument
  (:issue:`4102`, :issue:`4014`) in ``*.hist`` plotting methods
- Fixed bug in ``PeriodIndex.map`` where using ``str`` would return the str
  representation of the index (:issue:`4136`)
- Fixed test failure ``test_time_series_plot_color_with_empty_kwargs`` when
  using custom matplotlib default colors (:issue:`4345`)
- Fix running of stata IO tests. Now uses temporary files to write
  (:issue:`4353`)
- Fixed an issue where ``DataFrame.sum`` was slower than ``DataFrame.mean``
  for integer valued frames (:issue:`4365`)
- ``read_html`` tests now work with Python 2.6 (:issue:`4351`)
- Fixed bug where ``network`` testing was throwing ``NameError`` because a
  local variable was undefined (:issue:`4381`)
- In ``to_json``, raise if a passed ``orient`` would cause loss of data
  because of a duplicate index (:issue:`4359`)
- In ``to_json``, fix date handling so milliseconds are the default timestamp
  as the docstring says (:issue:`4362`).
- ``as_index`` is no longer ignored when doing groupby apply (:issue:`4648`,
  :issue:`3417`)
- JSON NaT handling fixed, NaTs are now serialized to `null` (:issue:`4498`)
- Fixed JSON handling of escapable characters in JSON object keys
  (:issue:`4593`)
- Fixed passing ``keep_default_na=False`` when ``na_values=None``
  (:issue:`4318`)
- Fixed bug with ``values`` raising an error on a DataFrame with duplicate
  columns and mixed dtypes, surfaced in (:issue:`4377`)
- Fixed bug with duplicate columns and type conversion in ``read_json`` when
  ``orient='split'`` (:issue:`4377`)
- Fixed JSON bug where locales with decimal separators other than '.' threw
  exceptions when encoding / decoding certain values. (:issue:`4918`)
- Fix ``.iat`` indexing with a ``PeriodIndex`` (:issue:`4390`)
- Fixed an issue where ``PeriodIndex`` joining with self was returning a new
  instance rather than the same instance (:issue:`4379`); also adds a test
  for this for the other index types
- Fixed a bug with all the dtypes being converted to object when using the
  CSV cparser with the usecols parameter (:issue:`3192`)
- Fix an issue in merging blocks where the resulting DataFrame had partially
  set _ref_locs (:issue:`4403`)
- Fixed an issue where hist subplots were being overwritten when they were
  called using the top level matplotlib API (:issue:`4408`)
- Fixed a bug where calling ``Series.astype(str)`` would truncate the string
  (:issue:`4405`, :issue:`4437`)
- Fixed a py3 compat issue where bytes were being repr'd as tuples
  (:issue:`4455`)
- Fixed Panel attribute naming conflict if item is named 'a'
  (:issue:`3440`)
- Fixed an issue where duplicate indexes were raising when plotting
  (:issue:`4486`)
- Fixed an issue where cumsum and cumprod didn't work with bool dtypes
  (:issue:`4170`, :issue:`4440`)
- Fixed Panel slicing issued in ``xs`` that was returning an incorrect dimmed
  object (:issue:`4016`)
- Fix resampling bug where custom reduce function not used if only one group
  (:issue:`3849`, :issue:`4494`)
- Fixed Panel assignment with a transposed frame (:issue:`3830`)
- Raise on set indexing with a Panel and a Panel as a value which needs
  alignment (:issue:`3777`)
- frozenset objects now raise in the ``Series`` constructor (:issue:`4482`,
  :issue:`4480`)
- Fixed issue with sorting a duplicate MultiIndex that has multiple dtypes
  (:issue:`4516`)
- Fixed bug in ``DataFrame.set_values`` which was causing name attributes to
  be lost when expanding the index. (:issue:`3742`, :issue:`4039`)
- Fixed issue where individual ``names``, ``levels`` and ``labels`` could be
  set on ``MultiIndex`` without validation (:issue:`3714`, :issue:`4039`)
- Fixed (:issue:`3334`) in pivot_table. Margins did not compute if values is
  the index.
- Fix bug in having a rhs of ``np.timedelta64`` or ``np.offsets.DateOffset``
  when operating with datetimes (:issue:`4532`)
- Fix arithmetic with series/datetimeindex and ``np.timedelta64`` not working
  the same (:issue:`4134`) and buggy timedelta in NumPy 1.6 (:issue:`4135`)
- Fix bug in ``pd.read_clipboard`` on windows with PY3 (:issue:`4561`); not
  decoding properly
- ``tslib.get_period_field()`` and ``tslib.get_period_field_arr()`` now raise
  if code argument out of range (:issue:`4519`, :issue:`4520`)
- Fix boolean indexing on an empty series loses index names (:issue:`4235`),
  infer_dtype works with empty arrays.
- Fix reindexing with multiple axes; if an axes match was not replacing the
  current axes, leading to a possible lazy frequency inference issue
  (:issue:`3317`)
- Fixed issue where ``DataFrame.apply`` was reraising exceptions incorrectly
  (causing the original stack trace to be truncated).
- Fix selection with ``ix/loc`` and non_unique selectors (:issue:`4619`)
- Fix assignment with iloc/loc involving a dtype change in an existing column
  (:issue:`4312`, :issue:`5702`) have internal setitem_with_indexer in core/indexing
  to use Block.setitem
- Fixed bug where thousands operator was not handled correctly for floating
  point numbers in csv_import (:issue:`4322`)
- Fix an issue with CacheableOffset not properly being used by many
  DateOffset; this prevented the DateOffset from being cached (:issue:`4609`)
- Fix boolean comparison with a DataFrame on the lhs, and a list/tuple on the
  rhs (:issue:`4576`)
- Fix error/dtype conversion with setitem of ``None`` on ``Series/DataFrame``
  (:issue:`4667`)
- Fix decoding based on a passed in non-default encoding in ``pd.read_stata``
  (:issue:`4626`)
- Fix ``DataFrame.from_records`` with a plain-vanilla ``ndarray``.
  (:issue:`4727`)
- Fix some inconsistencies with ``Index.rename`` and ``MultiIndex.rename``,
  etc. (:issue:`4718`, :issue:`4628`)
- Bug in using ``iloc/loc`` with a cross-sectional and duplicate indices
  (:issue:`4726`)
- Bug with using ``QUOTE_NONE`` with ``to_csv`` causing ``Exception``.
  (:issue:`4328`)
- Bug with Series indexing not raising an error when the right-hand-side has
  an incorrect length (:issue:`2702`)
- Bug in MultiIndexing with a partial string selection as one part of a
  MultIndex (:issue:`4758`)
- Bug with reindexing on the index with a non-unique index will now raise
  ``ValueError`` (:issue:`4746`)
- Bug in setting with ``loc/ix`` a single indexer with a MultiIndex axis and
  a NumPy array, related to (:issue:`3777`)
- Bug in concatenation with duplicate columns across dtypes not merging with
  axis=0 (:issue:`4771`, :issue:`4975`)
- Bug in ``iloc`` with a slice index failing (:issue:`4771`)
- Incorrect error message with no colspecs or width in ``read_fwf``.
  (:issue:`4774`)
- Fix bugs in indexing in a Series with a duplicate index (:issue:`4548`,
  :issue:`4550`)
- Fixed bug with reading compressed files with ``read_fwf`` in Python 3.
  (:issue:`3963`)
- Fixed an issue with a duplicate index and assignment with a dtype change
  (:issue:`4686`)
- Fixed bug with reading compressed files in as ``bytes`` rather than ``str``
  in Python 3. Simplifies bytes-producing file-handling in Python 3
  (:issue:`3963`, :issue:`4785`).
- Fixed an issue related to ticklocs/ticklabels with log scale bar plots
  across different versions of matplotlib (:issue:`4789`)
- Suppressed DeprecationWarning associated with internal calls issued by
  repr() (:issue:`4391`)
- Fixed an issue with a duplicate index and duplicate selector with ``.loc``
  (:issue:`4825`)
- Fixed an issue with ``DataFrame.sort_index`` where, when sorting by a
  single column and passing a list for ``ascending``, the argument for
  ``ascending`` was being interpreted as ``True`` (:issue:`4839`,
  :issue:`4846`)
- Fixed ``Panel.tshift`` not working. Added `freq` support to ``Panel.shift``
  (:issue:`4853`)
- Fix an issue in TextFileReader w/ Python engine (i.e. PythonParser)
  with thousands != "," (:issue:`4596`)
- Bug in getitem with a duplicate index when using where (:issue:`4879`)
- Fix Type inference code coerces float column into datetime (:issue:`4601`)
- Fixed ``_ensure_numeric`` does not check for complex numbers
  (:issue:`4902`)
- Fixed a bug in ``Series.hist`` where two figures were being created when
  the ``by`` argument was passed (:issue:`4112`, :issue:`4113`).
- Fixed a bug in ``convert_objects`` for > 2 ndims (:issue:`4937`)
- Fixed a bug in DataFrame/Panel cache insertion and subsequent indexing
  (:issue:`4939`, :issue:`5424`)
- Fixed string methods for ``FrozenNDArray`` and ``FrozenList``
  (:issue:`4929`)
- Fixed a bug with setting invalid or out-of-range values in indexing
  enlargement scenarios (:issue:`4940`)
- Tests for fillna on empty Series (:issue:`4346`), thanks @immerrr
- Fixed ``copy()`` to shallow copy axes/indices as well and thereby keep
  separate metadata. (:issue:`4202`, :issue:`4830`)
- Fixed skiprows option in Python parser for read_csv (:issue:`4382`)
- Fixed bug preventing ``cut`` from working with ``np.inf`` levels without
  explicitly passing labels (:issue:`3415`)
- Fixed wrong check for overlapping in ``DatetimeIndex.union``
  (:issue:`4564`)
- Fixed conflict between thousands separator and date parser in csv_parser
  (:issue:`4678`)
- Fix appending when dtypes are not the same (error showing mixing
  float/np.datetime64) (:issue:`4993`)
- Fix repr for DateOffset. No longer show duplicate entries in kwds.
  Removed unused offset fields. (:issue:`4638`)
- Fixed wrong index name during read_csv if using usecols. Applies to c
  parser only. (:issue:`4201`)
- ``Timestamp`` objects can now appear in the left hand side of a comparison
  operation with a ``Series`` or ``DataFrame`` object (:issue:`4982`).
- Fix a bug when indexing with ``np.nan`` via ``iloc/loc`` (:issue:`5016`)
- Fixed a bug where low memory c parser could create different types in
  different chunks of the same file. Now coerces to numerical type or raises
  warning. (:issue:`3866`)
- Fix a bug where reshaping a ``Series`` to its own shape raised
  ``TypeError`` (:issue:`4554`) and other reshaping issues.
- Bug in setting with ``ix/loc`` and a mixed int/string index (:issue:`4544`)
- Make sure series-series boolean comparisons are label based (:issue:`4947`)
- Bug in multi-level indexing with a Timestamp partial indexer
  (:issue:`4294`)
- Tests/fix for MultiIndex construction of an all-nan frame (:issue:`4078`)
- Fixed a bug where :func:`~pandas.read_html` wasn't correctly inferring
  values of tables with commas (:issue:`5029`)
- Fixed a bug where :func:`~pandas.read_html` wasn't providing a stable
  ordering of returned tables (:issue:`4770`, :issue:`5029`).
- Fixed a bug where :func:`~pandas.read_html` was incorrectly parsing when
  passed ``index_col=0`` (:issue:`5066`).
- Fixed a bug where :func:`~pandas.read_html` was incorrectly inferring the
  type of headers (:issue:`5048`).
- Fixed a bug where ``DatetimeIndex`` joins with ``PeriodIndex`` caused a
  stack overflow (:issue:`3899`).
- Fixed a bug where ``groupby`` objects didn't allow plots (:issue:`5102`).
- Fixed a bug where ``groupby`` objects weren't tab-completing column names
  (:issue:`5102`).
- Fixed a bug where ``groupby.plot()`` and friends were duplicating figures
  multiple times (:issue:`5102`).
- Provide automatic conversion of ``object`` dtypes on fillna, related
  (:issue:`5103`)
- Fixed a bug where default options were being overwritten in the option
  parser cleaning (:issue:`5121`).
- Treat a list/ndarray identically for ``iloc`` indexing with list-like
  (:issue:`5006`)
- Fix ``MultiIndex.get_level_values()`` with missing values (:issue:`5074`)
- Fix bound checking for Timestamp() with datetime64 input (:issue:`4065`)
- Fix a bug where ``TestReadHtml`` wasn't calling the correct ``read_html()``
  function (:issue:`5150`).
- Fix a bug with ``NDFrame.replace()`` which made replacement appear as
  though it was (incorrectly) using regular expressions (:issue:`5143`).
- Fix better error message for to_datetime (:issue:`4928`)
- Made sure different locales are tested on travis-ci (:issue:`4918`). Also
  adds a couple of utilities for getting locales and setting locales with a
  context manager.
- Fixed segfault on ``isnull(MultiIndex)`` (now raises an error instead)
  (:issue:`5123`, :issue:`5125`)
- Allow duplicate indices when performing operations that align
  (:issue:`5185`, :issue:`5639`)
- Compound dtypes in a constructor raise ``NotImplementedError``
  (:issue:`5191`)
- Bug in comparing duplicate frames (:issue:`4421`) related
- Bug in describe on duplicate frames
- Bug in ``to_datetime`` with a format and ``coerce=True`` not raising
  (:issue:`5195`)
- Bug in ``loc`` setting with multiple indexers and a rhs of a Series that
  needs broadcasting (:issue:`5206`)
- Fixed bug where inplace setting of levels or labels on ``MultiIndex`` would
  not clear cached ``values`` property and therefore return wrong ``values``.
  (:issue:`5215`)
- Fixed bug where filtering a grouped DataFrame or Series did not maintain
  the original ordering (:issue:`4621`).
- Fixed ``Period`` with a business date freq to always roll-forward if on a
  non-business date. (:issue:`5203`)
- Fixed bug in Excel writers where frames with duplicate column names weren't
  written correctly. (:issue:`5235`)
- Fixed issue with ``drop`` and a non-unique index on Series (:issue:`5248`)
- Fixed segfault in C parser caused by passing more names than columns in
  the file. (:issue:`5156`)
- Fix ``Series.isin`` with date/time-like dtypes (:issue:`5021`)
- C and Python Parser can now handle the more common MultiIndex column
  format which doesn't have a row for index names (:issue:`4702`)
- Bug when trying to use an out-of-bounds date as an object dtype
  (:issue:`5312`)
- Bug when trying to display an embedded PandasObject (:issue:`5324`)
- Allows operating of Timestamps to return a datetime if the result is out-of-bounds
  related (:issue:`5312`)
- Fix return value/type signature of ``initObjToJSON()`` to be compatible
  with numpy's ``import_array()`` (:issue:`5334`, :issue:`5326`)
- Bug when renaming then set_index on a DataFrame (:issue:`5344`)
- Test suite no longer leaves around temporary files when testing graphics. (:issue:`5347`)
  (thanks for catching this @yarikoptic!)
- Fixed html tests on win32. (:issue:`4580`)
- Make sure that ``head/tail`` are ``iloc`` based, (:issue:`5370`)
- Fixed bug for ``PeriodIndex`` string representation if there are 1 or 2
  elements. (:issue:`5372`)
- The GroupBy methods ``transform`` and ``filter`` can be used on Series
  and DataFrames that have repeated (non-unique) indices. (:issue:`4620`)
- Fix empty series not printing name in repr (:issue:`4651`)
- Make tests create temp files in temp directory by default. (:issue:`5419`)
- ``pd.to_timedelta`` of a scalar returns a scalar (:issue:`5410`)
- ``pd.to_timedelta`` accepts ``NaN`` and ``NaT``, returning ``NaT`` instead of raising (:issue:`5437`)
- performance improvements in ``isnull`` on larger size pandas objects
- Fixed various setitem with 1d ndarray that does not have a matching
  length to the indexer (:issue:`5508`)
- Bug in getitem with a MultiIndex and ``iloc`` (:issue:`5528`)
- Bug in delitem on a Series (:issue:`5542`)
- Bug fix in apply when using custom function and objects are not mutated (:issue:`5545`)
- Bug in selecting from a non-unique index with ``loc`` (:issue:`5553`)
- Bug in groupby returning non-consistent types when user function returns a ``None``, (:issue:`5592`)
- Work around regression in numpy 1.7.0 which erroneously raises IndexError from ``ndarray.item`` (:issue:`5666`)
- Bug in repeated indexing of object with resultant non-unique index (:issue:`5678`)
- Bug in fillna with Series and a passed series/dict (:issue:`5703`)
- Bug in groupby transform with a datetime-like grouper (:issue:`5712`)
- Bug in MultiIndex selection in PY3 when using certain keys (:issue:`5725`)
- Row-wise concat of differing dtypes failing in certain cases (:issue:`5754`)

pandas 0.12.0
-------------

**Release date:** 2013-07-24

New Features
~~~~~~~~~~~~

- ``pd.read_html()`` can now parse HTML strings, files or urls and returns a
  list of ``DataFrame`` s courtesy of @cpcloud. (:issue:`3477`,
  :issue:`3605`, :issue:`3606`)
- Support for reading Amazon S3 files. (:issue:`3504`)
- Added module for reading and writing JSON strings/files: pandas.io.json
  includes ``to_json`` DataFrame/Series method, and a ``read_json`` top-level reader
  various issues (:issue:`1226`, :issue:`3804`, :issue:`3876`, :issue:`3867`, :issue:`1305`)
- Added module for reading and writing Stata files: pandas.io.stata (:issue:`1512`)
  includes ``to_stata`` DataFrame method, and a ``read_stata`` top-level reader
- Added support for writing in ``to_csv`` and reading in ``read_csv``,
  MultiIndex columns. The ``header`` option in ``read_csv`` now accepts a
  list of the rows from which to read the index. Added the option,
  ``tupleize_cols`` to provide compatibility for the pre 0.12 behavior of
  writing and reading MultiIndex columns via a list of tuples. The default in
  0.12 is to write lists of tuples and *not* interpret list of tuples as a
  MultiIndex column.
  Note: The default value will change in 0.12 to make the default *to* write and
  read MultiIndex columns in the new format. (:issue:`3571`, :issue:`1651`, :issue:`3141`)
- Add iterator to ``Series.str`` (:issue:`3638`)
- ``pd.set_option()`` now allows N option, value pairs (:issue:`3667`).
- Added keyword parameters for different types of scatter_matrix subplots
- A ``filter`` method on grouped Series or DataFrames returns a subset of
  the original (:issue:`3680`, :issue:`919`)
- Access to historical Google Finance data in pandas.io.data (:issue:`3814`)
- DataFrame plotting methods can sample column colors from a Matplotlib
  colormap via the ``colormap`` keyword. (:issue:`3860`)

Improvements to existing features
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- Fixed various issues with internal pprinting code, the repr() for various objects
  including TimeStamp and Index now produces valid Python code strings and
  can be used to recreate the object, (:issue:`3038`, :issue:`3379`, :issue:`3251`, :issue:`3460`)
- ``convert_objects`` now accepts a ``copy`` parameter (defaults to ``True``)
- ``HDFStore``

  - will retain index attributes (freq,tz,name) on recreation (:issue:`3499`,:issue:`4098`)
  - will warn with a ``AttributeConflictWarning`` if you are attempting to append
    an index with a different frequency than the existing, or attempting
    to append an index with a different name than the existing
  - support datelike columns with a timezone as data_columns (:issue:`2852`)
  - table writing performance improvements.
  - support python3 (via ``PyTables 3.0.0``) (:issue:`3750`)
- Add modulo operator to Series, DataFrame
- Add ``date`` method to DatetimeIndex
- Add ``dropna`` argument to pivot_table (:issue: `3820`)
- Simplified the API and added a describe method to Categorical
- ``melt`` now accepts the optional parameters ``var_name`` and ``value_name``
  to specify custom column names of the returned DataFrame (:issue:`3649`),
  thanks @hoechenberger. If ``var_name`` is not specified and ``dataframe.columns.name``
  is not None, then this will be used as the ``var_name`` (:issue:`4144`).
  Also support for MultiIndex columns.
- clipboard functions use pyperclip (no dependencies on Windows, alternative
  dependencies offered for Linux) (:issue:`3837`).
- Plotting functions now raise a ``TypeError`` before trying to plot anything
  if the associated objects have a dtype of ``object`` (:issue:`1818`,
  :issue:`3572`, :issue:`3911`, :issue:`3912`), but they will try to convert object
  arrays to numeric arrays if possible so that you can still plot, for example, an
  object array with floats. This happens before any drawing takes place which
  eliminates any spurious plots from showing up.
- Added Faq section on repr display options, to help users customize their setup.
- ``where`` operations that result in block splitting are much faster (:issue:`3733`)
- Series and DataFrame hist methods now take a ``figsize`` argument (:issue:`3834`)
- DatetimeIndexes no longer try to convert mixed-integer indexes during join
  operations (:issue:`3877`)
- Add ``unit`` keyword to ``Timestamp`` and ``to_datetime`` to enable passing of
  integers or floats that are in an epoch unit of ``D, s, ms, us, ns``, thanks @mtkini (:issue:`3969`)
  (e.g. unix timestamps or epoch ``s``, with fractional seconds allowed) (:issue:`3540`)
- DataFrame corr method (spearman) is now cythonized.
- Improved ``network`` test decorator to catch ``IOError`` (and therefore
  ``URLError`` as well). Added ``with_connectivity_check`` decorator to allow
  explicitly checking a website as a proxy for seeing if there is network
  connectivity. Plus, new ``optional_args`` decorator factory for decorators.
  (:issue:`3910`, :issue:`3914`)
- ``read_csv`` will now throw a more informative error message when a file
  contains no columns, e.g., all newline characters
- Added ``layout`` keyword to DataFrame.hist() for more customizable layout (:issue:`4050`)
- Timestamp.min and Timestamp.max now represent valid Timestamp instances instead
  of the default datetime.min and datetime.max (respectively), thanks @SleepingPills
- ``read_html`` now raises when no tables are found and BeautifulSoup==4.2.0
  is detected (:issue:`4214`)

API Changes
~~~~~~~~~~~

- ``HDFStore``

  - When removing an object, ``remove(key)`` raises
    ``KeyError`` if the key is not a valid store object.
  - raise a ``TypeError`` on passing ``where`` or ``columns``
    to select with a Storer; these are invalid parameters at this time (:issue:`4189`)
  - can now specify an ``encoding`` option to ``append/put``
    to enable alternate encodings (:issue:`3750`)
  - enable support for ``iterator/chunksize`` with ``read_hdf``
- The repr() for (Multi)Index now obeys display.max_seq_items rather
  then NumPy threshold print options. (:issue:`3426`, :issue:`3466`)
- Added mangle_dupe_cols option to read_table/csv, allowing users
  to control legacy behaviour re dupe cols (A, A.1, A.2 vs A, A ) (:issue:`3468`)
  Note: The default value will change in 0.12 to the "no mangle" behaviour,
  If your code relies on this behaviour, explicitly specify mangle_dupe_cols=True
  in your calls.
- Do not allow astypes on ``datetime64[ns]`` except to ``object``, and
  ``timedelta64[ns]`` to ``object/int`` (:issue:`3425`)
- The behavior of ``datetime64`` dtypes has changed with respect to certain
  so-called reduction operations (:issue:`3726`). The following operations now
  raise a ``TypeError`` when performed on a ``Series`` and return an *empty*
  ``Series`` when performed on a ``DataFrame`` similar to performing these
  operations on, for example, a ``DataFrame`` of ``slice`` objects:
  - sum, prod, mean, std, var, skew, kurt, corr, and cov
- Do not allow datetimelike/timedeltalike creation except with valid types
  (e.g. cannot pass ``datetime64[ms]``) (:issue:`3423`)
- Add ``squeeze`` keyword to ``groupby`` to allow reduction from
  DataFrame -> Series if groups are unique. Regression from 0.10.1,
  partial revert on (:issue:`2893`) with (:issue:`3596`)
- Raise on ``iloc`` when boolean indexing with a label based indexer mask
  e.g. a boolean Series, even with integer labels, will raise. Since ``iloc``
  is purely positional based, the labels on the Series are not alignable (:issue:`3631`)
- The ``raise_on_error`` option to plotting methods is obviated by :issue:`3572`,
  so it is removed. Plots now always raise when data cannot be plotted or the
  object being plotted has a dtype of ``object``.
- ``DataFrame.interpolate()`` is now deprecated. Please use
  ``DataFrame.fillna()`` and ``DataFrame.replace()`` instead (:issue:`3582`,
  :issue:`3675`, :issue:`3676`).
- the ``method`` and ``axis`` arguments of ``DataFrame.replace()`` are
  deprecated
- ``DataFrame.replace`` 's ``infer_types`` parameter is removed and now
  performs conversion by default. (:issue:`3907`)
- Deprecated display.height, display.width is now only a formatting option
  does not control triggering of summary, similar to < 0.11.0.
- Add the keyword ``allow_duplicates`` to ``DataFrame.insert`` to allow a duplicate column
  to be inserted if ``True``, default is ``False`` (same as prior to 0.12) (:issue:`3679`)
- io API changes

  - added ``pandas.io.api`` for i/o imports
  - removed ``Excel`` support to ``pandas.io.excel``
  - added top-level ``pd.read_sql`` and ``to_sql`` DataFrame methods
  - removed ``clipboard`` support to ``pandas.io.clipboard``
  - replace top-level and instance methods ``save`` and ``load`` with
    top-level ``read_pickle`` and ``to_pickle`` instance method, ``save`` and
    ``load`` will give deprecation warning.
- the ``method`` and ``axis`` arguments of ``DataFrame.replace()`` are
  deprecated
- set FutureWarning to require data_source, and to replace year/month with
  expiry date in pandas.io options. This is in preparation to add options
  data from Google (:issue:`3822`)
- the ``method`` and ``axis`` arguments of ``DataFrame.replace()`` are
  deprecated
- Implement ``__nonzero__`` for ``NDFrame`` objects (:issue:`3691`, :issue:`3696`)
- ``as_matrix`` with mixed signed and unsigned dtypes will result in 2 x the lcd of the unsigned
  as an int, maxing with ``int64``, to avoid precision issues (:issue:`3733`)
- ``na_values`` in a list provided to ``read_csv/read_excel`` will match string and numeric versions
  e.g. ``na_values=['99']`` will match 99 whether the column ends up being int, float, or string (:issue:`3611`)
- ``read_html`` now defaults to ``None`` when reading, and falls back on
  ``bs4`` + ``html5lib`` when lxml fails to parse. a list of parsers to try
  until success is also valid
- more consistency in the to_datetime return types (give string/array of string inputs) (:issue:`3888`)
- The internal ``pandas`` class hierarchy has changed (slightly). The
  previous ``PandasObject`` now is called ``PandasContainer`` and a new
  ``PandasObject`` has become the base class for ``PandasContainer`` as well
  as ``Index``, ``Categorical``, ``GroupBy``, ``SparseList``, and
  ``SparseArray`` (+ their base classes). Currently, ``PandasObject``
  provides string methods (from ``StringMixin``). (:issue:`4090`, :issue:`4092`)
- New ``StringMixin`` that, given a ``__unicode__`` method, gets Python 2 and
  Python 3 compatible string methods (``__str__``, ``__bytes__``, and
  ``__repr__``). Plus string safety throughout. Now employed in many places
  throughout the pandas library. (:issue:`4090`, :issue:`4092`)

Experimental Features
~~~~~~~~~~~~~~~~~~~~~

- Added experimental ``CustomBusinessDay`` class to support ``DateOffsets``
  with custom holiday calendars and custom weekmasks. (:issue:`2301`)

Bug Fixes
~~~~~~~~~

- Fixed an esoteric excel reading bug, xlrd>= 0.9.0 now required for excel
  support. Should provide python3 support (for reading) which has been
  lacking. (:issue:`3164`)
- Disallow Series constructor called with MultiIndex which caused segfault (:issue:`4187`)
- Allow unioning of date ranges sharing a timezone (:issue:`3491`)
- Fix to_csv issue when having a large number of rows and ``NaT`` in some
  columns (:issue:`3437`)
- ``.loc`` was not raising when passed an integer list (:issue:`3449`)
- Unordered time series selection was misbehaving when using label slicing (:issue:`3448`)
- Fix sorting in a frame with a list of columns which contains datetime64[ns] dtypes (:issue:`3461`)
- DataFrames fetched via FRED now handle '.' as a NaN. (:issue:`3469`)
- Fix regression in a DataFrame apply with axis=1, objects were not being converted back
  to base dtypes correctly (:issue:`3480`)
- Fix issue when storing uint dtypes in an HDFStore. (:issue:`3493`)
- Non-unique index support clarified (:issue:`3468`)

  - Addressed handling of dupe columns in df.to_csv new and old (:issue:`3454`, :issue:`3457`)
  - Fix assigning a new index to a duplicate index in a DataFrame would fail (:issue:`3468`)
  - Fix construction of a DataFrame with a duplicate index
  - ref_locs support to allow duplicative indices across dtypes,
    allows iget support to always find the index (even across dtypes) (:issue:`2194`)
  - applymap on a DataFrame with a non-unique index now works
    (removed warning) (:issue:`2786`), and fix (:issue:`3230`)
  - Fix to_csv to handle non-unique columns (:issue:`3495`)
  - Duplicate indexes with getitem will return items in the correct order (:issue:`3455`, :issue:`3457`)
    and handle missing elements like unique indices (:issue:`3561`)
  - Duplicate indexes with and empty DataFrame.from_records will return a correct frame (:issue:`3562`)
  - Concat to produce a non-unique columns when duplicates are across dtypes is fixed (:issue:`3602`)
  - Non-unique indexing with a slice via ``loc`` and friends fixed (:issue:`3659`)
  - Allow insert/delete to non-unique columns (:issue:`3679`)
  - Extend ``reindex`` to correctly deal with non-unique indices (:issue:`3679`)
  - ``DataFrame.itertuples()`` now works with frames with duplicate column
    names (:issue:`3873`)
  - Bug in non-unique indexing via ``iloc`` (:issue:`4017`); added ``takeable`` argument to
    ``reindex`` for location-based taking
  - Allow non-unique indexing in series via ``.ix/.loc`` and ``__getitem__`` (:issue:`4246`)
  - Fixed non-unique indexing memory allocation issue with ``.ix/.loc`` (:issue:`4280`)

- Fixed bug in groupby with empty series referencing a variable before assignment. (:issue:`3510`)
- Allow index name to be used in groupby for non MultiIndex (:issue:`4014`)
- Fixed bug in mixed-frame assignment with aligned series (:issue:`3492`)
- Fixed bug in selecting month/quarter/year from a series would not select the time element
  on the last day (:issue:`3546`)
- Fixed a couple of MultiIndex rendering bugs in df.to_html() (:issue:`3547`, :issue:`3553`)
- Properly convert np.datetime64 objects in a Series (:issue:`3416`)
- Raise a ``TypeError`` on invalid datetime/timedelta operations
  e.g. add datetimes, multiple timedelta x datetime
- Fix ``.diff`` on datelike and timedelta operations (:issue:`3100`)
- ``combine_first`` not returning the same dtype in cases where it can (:issue:`3552`)
- Fixed bug with ``Panel.transpose`` argument aliases (:issue:`3556`)
- Fixed platform bug in ``PeriodIndex.take`` (:issue:`3579`)
- Fixed bud in incorrect conversion of datetime64[ns] in ``combine_first`` (:issue:`3593`)
- Fixed bug in reset_index with ``NaN`` in a MultiIndex (:issue:`3586`)
- ``fillna`` methods now raise a ``TypeError`` when the ``value`` parameter
  is a ``list`` or ``tuple``.
- Fixed bug where a time-series was being selected in preference to an actual column name
  in a frame (:issue:`3594`)
- Make secondary_y work properly for bar plots (:issue:`3598`)
- Fix modulo and integer division on Series,DataFrames to act similarly to ``float`` dtypes to return
  ``np.nan`` or ``np.inf`` as appropriate (:issue:`3590`)
- Fix incorrect dtype on groupby with ``as_index=False`` (:issue:`3610`)
- Fix ``read_csv/read_excel`` to correctly encode identical na_values, e.g. ``na_values=[-999.0,-999]``
  was failing (:issue:`3611`)
- Disable HTML output in qtconsole again. (:issue:`3657`)
- Reworked the new repr display logic, which users found confusing. (:issue:`3663`)
- Fix indexing issue in ndim >= 3 with ``iloc`` (:issue:`3617`)
- Correctly parse date columns with embedded (nan/NaT) into datetime64[ns] dtype in ``read_csv``
  when ``parse_dates`` is specified (:issue:`3062`)
- Fix not consolidating before to_csv (:issue:`3624`)
- Fix alignment issue when setitem in a DataFrame with a piece of a DataFrame (:issue:`3626`) or
  a mixed DataFrame and a Series (:issue:`3668`)
- Fix plotting of unordered DatetimeIndex (:issue:`3601`)
- ``sql.write_frame`` failing when writing a single column to sqlite (:issue:`3628`),
  thanks to @stonebig
- Fix pivoting with ``nan`` in the index (:issue:`3558`)
- Fix running of bs4 tests when it is not installed (:issue:`3605`)
- Fix parsing of html table (:issue:`3606`)
- ``read_html()`` now only allows a single backend: ``html5lib`` (:issue:`3616`)
- ``convert_objects`` with ``convert_dates='coerce'`` was parsing some single-letter strings into today's date
- ``DataFrame.from_records`` did not accept empty recarrays (:issue:`3682`)
- ``DataFrame.to_csv`` will succeed with the deprecated option ``nanRep``, @tdsmith
- ``DataFrame.to_html`` and ``DataFrame.to_latex`` now accept a path for
  their first argument (:issue:`3702`)
- Fix file tokenization error with \r delimiter and quoted fields (:issue:`3453`)
- Groupby transform with item-by-item not upcasting correctly (:issue:`3740`)
- Incorrectly read a HDFStore MultiIndex Frame with a column specification (:issue:`3748`)
- ``read_html`` now correctly skips tests (:issue:`3741`)
- PandasObjects raise TypeError when trying to hash (:issue:`3882`)
- Fix incorrect arguments passed to concat that are not list-like (e.g. concat(df1,df2)) (:issue:`3481`)
- Correctly parse when passed the ``dtype=str`` (or other variable-len string dtypes)
  in ``read_csv`` (:issue:`3795`)
- Fix index name not propagating when using ``loc/ix`` (:issue:`3880`)
- Fix groupby when applying a custom function resulting in a returned DataFrame was
  not converting dtypes (:issue:`3911`)
- Fixed a bug where ``DataFrame.replace`` with a compiled regular expression
  in the ``to_replace`` argument wasn't working (:issue:`3907`)
- Fixed ``__truediv__`` in Python 2.7 with ``numexpr`` installed to actually do true division when dividing
  two integer arrays with at least 10000 cells total (:issue:`3764`)
- Indexing with a string with seconds resolution not selecting from a time index (:issue:`3925`)
- csv parsers would loop infinitely if ``iterator=True`` but no ``chunksize`` was
  specified (:issue:`3967`), Python parser failing with ``chunksize=1``
- Fix index name not propagating when using ``shift``
- Fixed dropna=False being ignored with MultiIndex stack (:issue:`3997`)
- Fixed flattening of columns when renaming MultiIndex columns DataFrame (:issue:`4004`)
- Fix ``Series.clip`` for datetime series. NA/NaN threshold values will now throw ValueError (:issue:`3996`)
- Fixed insertion issue into DataFrame, after rename (:issue:`4032`)
- Fixed testing issue where too many sockets where open thus leading to a
  connection reset issue (:issue:`3982`, :issue:`3985`, :issue:`4028`,
  :issue:`4054`)
- Fixed failing tests in test_yahoo, test_google where symbols were not
  retrieved but were being accessed (:issue:`3982`, :issue:`3985`,
  :issue:`4028`, :issue:`4054`)
- ``Series.hist`` will now take the figure from the current environment if
  one is not passed
- Fixed bug where a 1xN DataFrame would barf on a 1xN mask (:issue:`4071`)
- Fixed running of ``tox`` under python3 where the pickle import was getting
  rewritten in an incompatible way (:issue:`4062`, :issue:`4063`)
- Fixed bug where sharex and sharey were not being passed to grouped_hist
  (:issue:`4089`)
- Fix bug where ``HDFStore`` will fail to append because of a different block
  ordering on-disk (:issue:`4096`)
- Better error messages on inserting incompatible columns to a frame (:issue:`4107`)
- Fixed bug in ``DataFrame.replace`` where a nested dict wasn't being
  iterated over when regex=False (:issue:`4115`)
- Fixed bug in ``convert_objects(convert_numeric=True)`` where a mixed numeric and
  object Series/Frame was not converting properly (:issue:`4119`)
- Fixed bugs in MultiIndex selection with column MultiIndex and duplicates
  (:issue:`4145`, :issue:`4146`)
- Fixed bug in the parsing of microseconds when using the ``format``
  argument in ``to_datetime`` (:issue:`4152`)
- Fixed bug in ``PandasAutoDateLocator`` where ``invert_xaxis`` triggered
  incorrectly ``MilliSecondLocator``  (:issue:`3990`)
- Fixed bug in ``Series.where`` where broadcasting a single element input vector
  to the length of the series resulted in multiplying the value
  inside the input (:issue:`4192`)
- Fixed bug in plotting that wasn't raising on invalid colormap for
  matplotlib 1.1.1 (:issue:`4215`)
- Fixed the legend displaying in ``DataFrame.plot(kind='kde')`` (:issue:`4216`)
- Fixed bug where Index slices weren't carrying the name attribute
  (:issue:`4226`)
- Fixed bug in initializing ``DatetimeIndex`` with an array of strings
  in a certain time zone (:issue:`4229`)
- Fixed bug where html5lib wasn't being properly skipped (:issue:`4265`)
- Fixed bug where get_data_famafrench wasn't using the correct file edges
  (:issue:`4281`)

pandas 0.11.0
-------------

**Release date:** 2013-04-22

New Features
~~~~~~~~~~~~

- New documentation section, ``10 Minutes to Pandas``
- New documentation section, ``Cookbook``
- Allow mixed dtypes (e.g ``float32/float64/int32/int16/int8``) to coexist in DataFrames and propagate in operations
- Add function to pandas.io.data for retrieving stock index components from Yahoo! finance (:issue:`2795`)
- Support slicing with time objects (:issue:`2681`)
- Added ``.iloc`` attribute, to support strict integer based indexing, analogous to ``.ix`` (:issue:`2922`)
- Added ``.loc`` attribute, to support strict label based indexing, analogous to ``.ix`` (:issue:`3053`)
- Added ``.iat`` attribute, to support fast scalar access via integers (replaces ``iget_value/iset_value``)
- Added ``.at`` attribute, to support fast scalar access via labels (replaces ``get_value/set_value``)
- Moved functionality from ``irow,icol,iget_value/iset_value`` to ``.iloc`` indexer (via ``_ixs`` methods in each object)
- Added support for expression evaluation using the ``numexpr`` library
- Added ``convert=boolean`` to ``take`` routines to translate negative indices to positive, defaults to True
- Added to_series() method to indices, to facilitate the creation of indexers (:issue:`3275`)

Improvements to existing features
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- Improved performance of df.to_csv() by up to 10x in some cases. (:issue:`3059`)
- added ``blocks`` attribute to DataFrames, to return a dict of dtypes to homogeneously dtyped DataFrames
- added keyword ``convert_numeric`` to ``convert_objects()`` to try to convert object dtypes to numeric types (default is False)
- ``convert_dates`` in ``convert_objects`` can now be ``coerce`` which will return
  a datetime64[ns] dtype with non-convertibles set as ``NaT``; will preserve an all-nan object
  (e.g. strings), default is True (to perform soft-conversion
- Series print output now includes the dtype by default
- Optimize internal reindexing routines (:issue:`2819`, :issue:`2867`)
- ``describe_option()`` now reports the default and current value of options.
- Add ``format`` option to ``pandas.to_datetime`` with faster conversion of strings that can be parsed with datetime.strptime
- Add ``axes`` property to ``Series`` for compatibility
- Add ``xs`` function to ``Series`` for compatibility
- Allow setitem in a frame where only mixed numerics are present (e.g. int and float), (:issue:`3037`)
- ``HDFStore``

  - Provide dotted attribute access to ``get`` from stores (e.g. store.df == store['df'])
  - New keywords ``iterator=boolean``, and ``chunksize=number_in_a_chunk`` are provided to support iteration on ``select`` and ``select_as_multiple`` (:issue:`3076`)
  - support ``read_hdf/to_hdf`` API similar to ``read_csv/to_csv`` (:issue:`3222`)

- Add ``squeeze`` method to possibly remove length 1 dimensions from an object.

  .. ipython:: python

     p = pd.Panel(np.random.randn(3,4,4),items=['ItemA','ItemB','ItemC'],
                  major_axis=pd.date_range('20010102',periods=4),
                  minor_axis=['A','B','C','D'])
     p
     p.reindex(items=['ItemA']).squeeze()
     p.reindex(items=['ItemA'],minor=['B']).squeeze()

- Improvement to Yahoo API access in ``pd.io.data.Options`` (:issue:`2758`)
- added option `display.max_seq_items` to control the number of elements printed per sequence pprinting it. (:issue:`2979`)
- added option `display.chop_threshold` to control display of small numerical values. (:issue:`2739`)
- added option `display.max_info_rows` to prevent verbose_info from being
  calculated for frames above 1M rows (configurable). (:issue:`2807`, :issue:`2918`)
- value_counts() now accepts a "normalize" argument, for normalized histograms. (:issue:`2710`).
- DataFrame.from_records now accepts not only dicts but any instance of the collections.Mapping ABC.
- Allow selection semantics via a string with a datelike index to work in both Series and DataFrames (:issue:`3070`)

  .. ipython:: python

      idx = pd.date_range("2001-10-1", periods=5, freq='M')
      ts = pd.Series(np.random.rand(len(idx)),index=idx)
      ts['2001']

      df = pd.DataFrame(dict(A = ts))
      df['2001']

- added option `display.mpl_style` providing a sleeker visual style for plots. Based on https://gist.github.com/huyng/816622 (:issue:`3075`).
- Improved performance across several core functions by taking memory ordering of
  arrays into account. Courtesy of @stephenwlin (:issue:`3130`)
- Improved performance of groupby transform method (:issue:`2121`)
- Handle "ragged" CSV files missing trailing delimiters in rows with missing fields
  when also providing explicit list of column names (so the parser knows how many columns to expect in the result) (:issue:`2981`)
- On a mixed DataFrame, allow setting with indexers with ndarray/DataFrame on rhs (:issue:`3216`)
- Treat boolean values as integers (values 1 and 0) for numeric operations. (:issue:`2641`)
- Add ``time`` method to DatetimeIndex (:issue:`3180`)
- Return NA when using Series.str[...] for values that are not long enough (:issue:`3223`)
- Display cursor coordinate information in time-series plots (:issue:`1670`)
- to_html() now accepts an optional "escape" argument to control reserved HTML character
  escaping (enabled by default) and escapes ``&``, in addition to ``<`` and ``>``.  (:issue:`2919`)

API Changes
~~~~~~~~~~~

- Do not automatically upcast numeric specified dtypes to ``int64`` or
  ``float64`` (:issue:`622` and :issue:`797`)
- DataFrame construction of lists and scalars, with no dtype present, will
  result in casting to ``int64`` or ``float64``, regardless of platform.
  This is not an apparent change in the API, but noting it.
- Guarantee that ``convert_objects()`` for Series/DataFrame always returns a
  copy
- groupby operations will respect dtypes for numeric float operations
  (float32/float64); other types will be operated on, and will try to cast
  back to the input dtype (e.g. if an int is passed, as long as the output
  doesn't have nans, then an int will be returned)
- backfill/pad/take/diff/ohlc will now support ``float32/int16/int8``
  operations
- Block types will upcast as needed in where/masking operations (:issue:`2793`)
- Series now automatically will try to set the correct dtype based on passed
  datetimelike objects (datetime/Timestamp)

  - timedelta64 are returned in appropriate cases (e.g. Series - Series,
    when both are datetime64)
  - mixed datetimes and objects (:issue:`2751`) in a constructor will be cast
    correctly
  - astype on datetimes to object are now handled (as well as NaT
    conversions to np.nan)
  - all timedelta like objects will be correctly assigned to ``timedelta64``
    with mixed ``NaN`` and/or ``NaT`` allowed

- arguments to DataFrame.clip were inconsistent to NumPy and Series clipping
  (:issue:`2747`)
- util.testing.assert_frame_equal now checks the column and index names (:issue:`2964`)
- Constructors will now return a more informative ValueError on failures
  when invalid shapes are passed
- Don't suppress TypeError in GroupBy.agg (:issue:`3238`)
- Methods return None when inplace=True (:issue:`1893`)
- ``HDFStore``

   - added the method ``select_column`` to select a single column from a table as a Series.
   - deprecated the ``unique`` method, can be replicated by ``select_column(key,column).unique()``
   - ``min_itemsize`` parameter will now automatically create data_columns for passed keys

- Downcast on pivot if possible (:issue:`3283`), adds argument ``downcast`` to ``fillna``
- Introduced options `display.height/width` for explicitly specifying terminal
  height/width in characters. Deprecated display.line_width, now replaced by display.width.
  These defaults are in effect for scripts as well, so unless disabled, previously
  very wide output will now be output as "expand_repr" style wrapped output.
- Various defaults for options (including display.max_rows) have been revised,
  after a brief survey concluded they were wrong for everyone. Now at w=80,h=60.
- HTML repr output in IPython qtconsole is once again controlled by the option
  `display.notebook_repr_html`, and on by default.

Bug Fixes
~~~~~~~~~

- Fix seg fault on empty data frame when fillna with ``pad`` or ``backfill``
  (:issue:`2778`)
- Single element ndarrays of datetimelike objects are handled
  (e.g. np.array(datetime(2001,1,1,0,0))), w/o dtype being passed
- 0-dim ndarrays with a passed dtype are handled correctly
  (e.g. np.array(0.,dtype='float32'))
- Fix some boolean indexing inconsistencies in Series.__getitem__/__setitem__
  (:issue:`2776`)
- Fix issues with DataFrame and Series constructor with integers that
  overflow ``int64`` and some mixed typed type lists (:issue:`2845`)

- ``HDFStore``

  - Fix weird PyTables error when using too many selectors in a where
    also correctly filter on any number of values in a Term expression
    (so not using numexpr filtering, but isin filtering)
  - Internally, change all variables to be private-like (now have leading
    underscore)
  - Fixes for query parsing to correctly interpret boolean and != (:issue:`2849`, :issue:`2973`)
  - Fixes for pathological case on SparseSeries with 0-len array and
    compression (:issue:`2931`)
  - Fixes bug with writing rows if part of a block was all-nan (:issue:`3012`)
  - Exceptions are now ValueError or TypeError as needed
  - A table will now raise if min_itemsize contains fields which are not queryables

- Bug showing up in applymap where some object type columns are converted (:issue:`2909`)
  had an incorrect default in convert_objects

- TimeDeltas

  - Series ops with a Timestamp on the rhs was throwing an exception (:issue:`2898`)
    added tests for Series ops with datetimes,timedeltas,Timestamps, and datelike
    Series on both lhs and rhs
  - Fixed subtle timedelta64 inference issue on py3 & NumPy 1.7.0 (:issue:`3094`)
  - Fixed some formatting issues on timedelta when negative
  - Support null checking on timedelta64, representing (and formatting) with NaT
  - Support setitem with np.nan value, converts to NaT
  - Support min/max ops in a Dataframe (abs not working, nor do we error on non-supported ops)
  - Support idxmin/idxmax/abs/max/min in a Series (:issue:`2989`, :issue:`2982`)

- Bug on in-place putmasking on an ``integer`` series that needs to be converted to
  ``float`` (:issue:`2746`)
- Bug in argsort of ``datetime64[ns]`` Series with ``NaT`` (:issue:`2967`)
- Bug in value_counts of ``datetime64[ns]`` Series (:issue:`3002`)
- Fixed printing of ``NaT`` in an index
- Bug in idxmin/idxmax of ``datetime64[ns]`` Series with ``NaT`` (:issue:`2982`)
- Bug in ``icol, take`` with negative indices was producing incorrect return
  values (see :issue:`2922`, :issue:`2892`), also check for out-of-bounds indices (:issue:`3029`)
- Bug in DataFrame column insertion when the column creation fails, existing frame is left in
  an irrecoverable state (:issue:`3010`)
- Bug in DataFrame update, combine_first where non-specified values could cause
  dtype changes (:issue:`3016`, :issue:`3041`)
- Bug in groupby with first/last where dtypes could change (:issue:`3041`, :issue:`2763`)
- Formatting of an index that has ``nan`` was inconsistent or wrong (would fill from
  other values), (:issue:`2850`)
- Unstack of a frame with no nans would always cause dtype upcasting (:issue:`2929`)
- Fix scalar datetime.datetime parsing bug in read_csv (:issue:`3071`)
- Fixed slow printing of large Dataframes, due to inefficient dtype
  reporting (:issue:`2807`)
- Fixed a segfault when using a function as grouper in groupby (:issue:`3035`)
- Fix pretty-printing of infinite data structures (closes :issue:`2978`)
- Fixed exception when plotting timeseries bearing a timezone (closes :issue:`2877`)
- str.contains ignored na argument (:issue:`2806`)
- Substitute warning for segfault when grouping with categorical grouper
  of mismatched length (:issue:`3011`)
- Fix exception in SparseSeries.density (:issue:`2083`)
- Fix upsampling bug with closed='left' and daily to daily data (:issue:`3020`)
- Fixed missing tick bars on scatter_matrix plot (:issue:`3063`)
- Fixed bug in Timestamp(d,tz=foo) when d is date() rather then datetime() (:issue:`2993`)
- series.plot(kind='bar') now respects pylab color scheme (:issue:`3115`)
- Fixed bug in reshape if not passed correct input, now raises TypeError (:issue:`2719`)
- Fixed a bug where Series ctor did not respect ordering if OrderedDict passed in (:issue:`3282`)
- Fix NameError issue on RESO_US (:issue:`2787`)
- Allow selection in an *unordered* timeseries to work similarly
  to an *ordered* timeseries (:issue:`2437`).
- Fix implemented ``.xs`` when called with ``axes=1`` and a level parameter (:issue:`2903`)
- Timestamp now supports the class method fromordinal similar to datetimes (:issue:`3042`)
- Fix issue with indexing a series with a boolean key and specifying a 1-len list on the rhs (:issue:`2745`)
  or a list on the rhs (:issue:`3235`)
- Fixed bug in groupby apply when kernel generate list of arrays having unequal len (:issue:`1738`)
- fixed handling of rolling_corr with center=True which could produce corr>1 (:issue:`3155`)
- Fixed issues where indices can be passed as 'index/column' in addition to 0/1 for the axis parameter
- PeriodIndex.tolist now boxes to Period (:issue:`3178`)
- PeriodIndex.get_loc KeyError now reports Period instead of ordinal (:issue:`3179`)
- df.to_records bug when handling MultiIndex (GH3189)
- Fix Series.__getitem__ segfault when index less than -length (:issue:`3168`)
- Fix bug when using Timestamp as a date parser (:issue:`2932`)
- Fix bug creating date range from Timestamp with time zone and passing same
  time zone (:issue:`2926`)
- Add comparison operators to Period object (:issue:`2781`)
- Fix bug when concatenating two Series into a DataFrame when they have the
  same name (:issue:`2797`)
- Fix automatic color cycling when plotting consecutive timeseries
  without color arguments (:issue:`2816`)
- fixed bug in the pickling of PeriodIndex (:issue:`2891`)
- Upcast/split blocks when needed in a mixed DataFrame when setitem
  with an indexer (:issue:`3216`)
- Invoking df.applymap on a dataframe with dupe cols now raises a ValueError (:issue:`2786`)
- Apply with invalid returned indices raise correct Exception (:issue:`2808`)
- Fixed a bug in plotting log-scale bar plots (:issue:`3247`)
- df.plot() grid on/off now obeys the mpl default style, just like
  series.plot(). (:issue:`3233`)
- Fixed a bug in the legend of plotting.andrews_curves() (:issue:`3278`)
- Produce a series on apply if we only generate a singular series and have
  a simple index (:issue:`2893`)
- Fix Python ASCII file parsing when integer falls outside of floating point
  spacing (:issue:`3258`)
- fixed pretty printing of sets (:issue:`3294`)
- Panel() and Panel.from_dict() now respects ordering when give OrderedDict (:issue:`3303`)
- DataFrame where with a datetimelike incorrectly selecting (:issue:`3311`)
- Ensure index casts work even in Int64Index
- Fix set_index segfault when passing MultiIndex (:issue:`3308`)
- Ensure pickles created in py2 can be read in py3
- Insert ellipsis in MultiIndex summary repr (:issue:`3348`)
- Groupby will handle mutation among an input groups columns (and fallback
  to non-fast apply) (:issue:`3380`)
- Eliminated unicode errors on FreeBSD when using MPL GTK backend (:issue:`3360`)
- Period.strftime should return unicode strings always (:issue:`3363`)
- Respect passed read_* chunksize in get_chunk function (:issue:`3406`)

pandas 0.10.1
-------------

**Release date:** 2013-01-22

New Features
~~~~~~~~~~~~

- Add data interface to World Bank WDI pandas.io.wb (:issue:`2592`)

API Changes
~~~~~~~~~~~

- Restored inplace=True behavior returning self (same object) with
  deprecation warning until 0.11 (:issue:`1893`)
- ``HDFStore``

  - refactored HFDStore to deal with non-table stores as objects, will allow future enhancements
  - removed keyword ``compression`` from ``put`` (replaced by keyword
    ``complib`` to be consistent across library)
  - warn `PerformanceWarning` if you are attempting to store types that will be pickled by PyTables

Improvements to existing features
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- ``HDFStore``

  - enables storing of MultiIndex dataframes (closes :issue:`1277`)
  - support data column indexing and selection, via ``data_columns`` keyword
    in append
  - support write chunking to reduce memory footprint, via ``chunksize``
    keyword to append
  - support automagic indexing via ``index`` keyword to append
  - support ``expectedrows`` keyword in append to inform ``PyTables`` about
    the expected table size
  - support ``start`` and ``stop`` keywords in select to limit the row
    selection space
  - added ``get_store`` context manager to automatically import with pandas
  - added column filtering via ``columns`` keyword in select
  - added methods append_to_multiple/select_as_multiple/select_as_coordinates
    to do multiple-table append/selection
  - added support for datetime64 in columns
  - added method ``unique`` to select the unique values in an indexable or
    data column
  - added method ``copy`` to copy an existing store (and possibly upgrade)
  - show the shape of the data on disk for non-table stores when printing the
    store
  - added ability to read PyTables flavor tables (allows compatibility to
    other HDF5 systems)

- Add ``logx`` option to DataFrame/Series.plot (:issue:`2327`, :issue:`2565`)
- Support reading gzipped data from file-like object
- ``pivot_table`` aggfunc can be anything used in GroupBy.aggregate (:issue:`2643`)
- Implement DataFrame merges in case where set cardinalities might overflow
  64-bit integer (:issue:`2690`)
- Raise exception in C file parser if integer dtype specified and have NA
  values. (:issue:`2631`)
- Attempt to parse ISO8601 format dates when parse_dates=True in read_csv for
  major performance boost in such cases (:issue:`2698`)
- Add methods ``neg`` and ``inv`` to Series
- Implement ``kind`` option in ``ExcelFile`` to indicate whether it's an XLS
  or XLSX file (:issue:`2613`)
- Documented a fast-path in pd.read_csv when parsing iso8601 datetime strings
  yielding as much as a 20x speedup.  (:issue:`5993`)


Bug Fixes
~~~~~~~~~

- Fix read_csv/read_table multithreading issues (:issue:`2608`)
- ``HDFStore``

  - correctly handle ``nan`` elements in string columns; serialize via the
    ``nan_rep`` keyword to append
  - raise correctly on non-implemented column types (unicode/date)
  - handle correctly ``Term`` passed types (e.g. ``index<1000``, when index
    is ``Int64``), (closes :issue:`512`)
  - handle Timestamp correctly in data_columns (closes :issue:`2637`)
  - contains correctly matches on non-natural names
  - correctly store ``float32`` dtypes in tables (if not other float types in
    the same table)

- Fix DataFrame.info bug with UTF8-encoded columns. (:issue:`2576`)
- Fix DatetimeIndex handling of FixedOffset tz (:issue:`2604`)
- More robust detection of being in IPython session for wide DataFrame
  console formatting (:issue:`2585`)
- Fix platform issues with ``file:///`` in unit test (:issue:`2564`)
- Fix bug and possible segfault when grouping by hierarchical level that
  contains NA values (:issue:`2616`)
- Ensure that MultiIndex tuples can be constructed with NAs (:issue:`2616`)
- Fix int64 overflow issue when unstacking MultiIndex with many levels
  (:issue:`2616`)
- Exclude non-numeric data from DataFrame.quantile by default (:issue:`2625`)
- Fix a Cython C int64 boxing issue causing read_csv to return incorrect
  results (:issue:`2599`)
- Fix groupby summing performance issue on boolean data (:issue:`2692`)
- Don't bork Series containing datetime64 values with to_datetime (:issue:`2699`)
- Fix DataFrame.from_records corner case when passed columns, index column,
  but empty record list (:issue:`2633`)
- Fix C parser-tokenizer bug with trailing fields. (:issue:`2668`)
- Don't exclude non-numeric data from GroupBy.max/min (:issue:`2700`)
- Don't lose time zone when calling DatetimeIndex.drop (:issue:`2621`)
- Fix setitem on a Series with a boolean key and a non-scalar as value
  (:issue:`2686`)
- Box datetime64 values in Series.apply/map (:issue:`2627`, :issue:`2689`)
- Up convert datetime + datetime64 values when concatenating frames (:issue:`2624`)
- Raise a more helpful error message in merge operations when one DataFrame
  has duplicate columns (:issue:`2649`)
- Fix partial date parsing issue occurring only when code is run at EOM
  (:issue:`2618`)
- Prevent MemoryError when using counting sort in sortlevel with
  high-cardinality MultiIndex objects (:issue:`2684`)
- Fix Period resampling bug when all values fall into a single bin (:issue:`2070`)
- Fix buggy interaction with usecols argument in read_csv when there is an
  implicit first index column (:issue:`2654`)
- Fix bug in ``Index.summary()`` where string format methods were being called incorrectly.
  (:issue:`3869`)

pandas 0.10.0
-------------

**Release date:** 2012-12-17

New Features
~~~~~~~~~~~~

- Brand new high-performance delimited file parsing engine written in C and
  Cython. 50% or better performance in many standard use cases with a
  fraction as much memory usage. (:issue:`407`, :issue:`821`)
- Many new file parser (read_csv, read_table) features:

  - Support for on-the-fly gzip or bz2 decompression (`compression` option)
  - Ability to get back numpy.recarray instead of DataFrame
    (`as_recarray=True`)
  - `dtype` option: explicit column dtypes
  - `usecols` option: specify list of columns to be read from a file. Good
    for reading very wide files with many irrelevant columns (:issue:`1216` :issue:`926`, :issue:`2465`)
  - Enhanced unicode decoding support via `encoding` option
  - `skipinitialspace` dialect option
  - Can specify strings to be recognized as True (`true_values`) or False
    (`false_values`)
  - High-performance `delim_whitespace` option for whitespace-delimited
    files; a preferred alternative to the '\s+' regular expression delimiter
  - Option to skip "bad" lines (wrong number of fields) that would otherwise
    have caused an error in the past (`error_bad_lines` and `warn_bad_lines`
    options)
  - Substantially improved performance in the parsing of integers with
    thousands markers and lines with comments
  - Easy of European (and other) decimal formats (`decimal` option) (:issue:`584`, :issue:`2466`)
  - Custom line terminators (e.g. lineterminator='~') (:issue:`2457`)
  - Handling of no trailing commas in CSV files (:issue:`2333`)
  - Ability to handle fractional seconds in date_converters (:issue:`2209`)
  - read_csv allow scalar arg to na_values (:issue:`1944`)
  - Explicit column dtype specification in read_* functions (:issue:`1858`)
  - Easier CSV dialect specification (:issue:`1743`)
  - Improve parser performance when handling special characters (:issue:`1204`)

- Google Analytics API integration with easy oauth2 workflow (:issue:`2283`)
- Add error handling to Series.str.encode/decode (:issue:`2276`)
- Add ``where`` and ``mask`` to Series (:issue:`2337`)
- Grouped histogram via `by` keyword in Series/DataFrame.hist (:issue:`2186`)
- Support optional ``min_periods`` keyword in ``corr`` and ``cov``
  for both Series and DataFrame (:issue:`2002`)
- Add ``duplicated`` and ``drop_duplicates`` functions to Series (:issue:`1923`)
- Add docs for ``HDFStore table`` format
- 'density' property in `SparseSeries` (:issue:`2384`)
- Add ``ffill`` and ``bfill`` convenience functions for forward- and
  backfilling time series data (:issue:`2284`)
- New option configuration system and functions `set_option`, `get_option`,
  `describe_option`, and `reset_option`. Deprecate `set_printoptions` and
  `reset_printoptions` (:issue:`2393`).
  You can also access options as attributes via ``pandas.options.X``
- Wide DataFrames can be viewed more easily in the console with new
  `expand_frame_repr` and `line_width` configuration options. This is on by
  default now (:issue:`2436`)
- Scikits.timeseries-like moving window functions via ``rolling_window`` (:issue:`1270`)

Experimental Features
~~~~~~~~~~~~~~~~~~~~~

- Add support for Panel4D, a named 4 Dimensional structure
- Add support for ndpanel factory functions, to create custom,
  domain-specific N-Dimensional containers

API Changes
~~~~~~~~~~~

- The default binning/labeling behavior for ``resample`` has been changed to
  `closed='left', label='left'` for daily and lower frequencies. This had
  been a large source of confusion for users. See "what's new" page for more
  on this. (:issue:`2410`)
- Methods with ``inplace`` option now return None instead of the calling
  (modified) object (:issue:`1893`)
- The special case DataFrame - TimeSeries doing column-by-column broadcasting
  has been deprecated. Users should explicitly do e.g. df.sub(ts, axis=0)
  instead. This is a legacy hack and can lead to subtle bugs.
- inf/-inf are no longer considered as NA by isnull/notnull. To be clear, this
  is legacy cruft from early pandas. This behavior can be globally re-enabled
  using the new option ``mode.use_inf_as_null`` (:issue:`2050`, :issue:`1919`)
- ``pandas.merge`` will now default to ``sort=False``. For many use cases
  sorting the join keys is not necessary, and doing it by default is wasteful
- Specify ``header=0`` explicitly to replace existing column names in file in
  read_* functions.
- Default column names for header-less parsed files (yielded by read_csv,
  etc.) are now the integers 0, 1, .... A new argument `prefix` has been
  added; to get the v0.9.x behavior specify ``prefix='X'`` (:issue:`2034`). This API
  change was made to make the default column names more consistent with the
  DataFrame constructor's default column names when none are specified.
- DataFrame selection using a boolean frame now preserves input shape
- If function passed to Series.apply yields a Series, result will be a
  DataFrame (:issue:`2316`)
- Values like YES/NO/yes/no will not be considered as boolean by default any
  longer in the file parsers. This can be customized using the new
  ``true_values`` and ``false_values`` options (:issue:`2360`)
- `obj.fillna()` is no longer valid; make `method='pad'` no longer the
  default option, to be more explicit about what kind of filling to
  perform. Add `ffill/bfill` convenience functions per above (:issue:`2284`)
- `HDFStore.keys()` now returns an absolute path-name for each key
- `to_string()` now always returns a unicode string. (:issue:`2224`)
- File parsers will not handle NA sentinel values arising from passed
  converter functions

Improvements to existing features
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- Add ``nrows`` option to DataFrame.from_records for iterators (:issue:`1794`)
- Unstack/reshape algorithm rewrite to avoid high memory use in cases where
  the number of observed key-tuples is much smaller than the total possible
  number that could occur (:issue:`2278`). Also improves performance in most cases.
- Support duplicate columns in DataFrame.from_records (:issue:`2179`)
- Add ``normalize`` option to Series/DataFrame.asfreq (:issue:`2137`)
- SparseSeries and SparseDataFrame construction from empty and scalar
  values now no longer create dense ndarrays unnecessarily (:issue:`2322`)
- ``HDFStore`` now supports hierarchical keys (:issue:`2397`)
- Support multiple query selection formats for ``HDFStore tables`` (:issue:`1996`)
- Support ``del store['df']`` syntax to delete HDFStores
- Add multi-dtype support for ``HDFStore tables``
- ``min_itemsize`` parameter can be specified in ``HDFStore table`` creation
- Indexing support in ``HDFStore tables`` (:issue:`698`)
- Add `line_terminator` option to DataFrame.to_csv (:issue:`2383`)
- added implementation of str(x)/unicode(x)/bytes(x) to major pandas data
  structures, which should do the right thing on both py2.x and py3.x. (:issue:`2224`)
- Reduce groupby.apply overhead substantially by low-level manipulation of
  internal NumPy arrays in DataFrames (:issue:`535`)
- Implement ``value_vars`` in ``melt`` and add ``melt`` to pandas namespace
  (:issue:`2412`)
- Added boolean comparison operators to Panel
- Enable ``Series.str.strip/lstrip/rstrip`` methods to take an argument (:issue:`2411`)
- The DataFrame ctor now respects column ordering when given
  an OrderedDict (:issue:`2455`)
- Assigning DatetimeIndex to Series changes the class to TimeSeries (:issue:`2139`)
- Improve performance of .value_counts method on non-integer data (:issue:`2480`)
- ``get_level_values`` method for MultiIndex return Index instead of ndarray (:issue:`2449`)
- ``convert_to_r_dataframe`` conversion for datetime values (:issue:`2351`)
- Allow ``DataFrame.to_csv`` to represent inf and nan differently (:issue:`2026`)
- Add ``min_i`` argument to ``nancorr`` to specify minimum required observations (:issue:`2002`)
- Add ``inplace`` option to ``sortlevel`` / ``sort`` functions on DataFrame (:issue:`1873`)
- Enable DataFrame to accept scalar constructor values like Series (:issue:`1856`)
- DataFrame.from_records now takes optional ``size`` parameter (:issue:`1794`)
- include iris dataset (:issue:`1709`)
- No datetime64 DataFrame column conversion of datetime.datetime with tzinfo (:issue:`1581`)
- Micro-optimizations in DataFrame for tracking state of internal consolidation (:issue:`217`)
- Format parameter in DataFrame.to_csv (:issue:`1525`)
- Partial string slicing for ``DatetimeIndex`` for daily and higher frequencies (:issue:`2306`)
- Implement ``col_space`` parameter in ``to_html`` and ``to_string`` in DataFrame (:issue:`1000`)
- Override ``Series.tolist`` and box datetime64 types (:issue:`2447`)
- Optimize ``unstack`` memory usage by compressing indices (:issue:`2278`)
- Fix HTML repr in IPython qtconsole if opening window is small (:issue:`2275`)
- Escape more special characters in console output (:issue:`2492`)
- df.select now invokes bool on the result of crit(x) (:issue:`2487`)

Bug Fixes
~~~~~~~~~

- Fix major performance regression in DataFrame.iteritems (:issue:`2273`)
- Fixes bug when negative period passed to Series/DataFrame.diff (:issue:`2266`)
- Escape tabs in console output to avoid alignment issues (:issue:`2038`)
- Properly box datetime64 values when retrieving cross-section from
  mixed-dtype DataFrame (:issue:`2272`)
- Fix concatenation bug leading to :issue:`2057`, :issue:`2257`
- Fix regression in Index console formatting (:issue:`2319`)
- Box Period data when assigning PeriodIndex to frame column (:issue:`2243`, :issue:`2281`)
- Raise exception on calling reset_index on Series with inplace=True (:issue:`2277`)
- Enable setting multiple columns in DataFrame with hierarchical columns
  (:issue:`2295`)
- Respect dtype=object in DataFrame constructor (:issue:`2291`)
- Fix DatetimeIndex.join bug with tz-aware indexes and how='outer' (:issue:`2317`)
- pop(...) and del works with DataFrame with duplicate columns (:issue:`2349`)
- Treat empty strings as NA in date parsing (rather than let dateutil do
  something weird) (:issue:`2263`)
- Prevent uint64 -> int64 overflows (:issue:`2355`)
- Enable joins between MultiIndex and regular Index (:issue:`2024`)
- Fix time zone metadata issue when unioning non-overlapping DatetimeIndex
  objects (:issue:`2367`)
- Raise/handle int64 overflows in parsers (:issue:`2247`)
- Deleting of consecutive rows in ``HDFStore tables``` is much faster than before
- Appending on a HDFStore would fail if the table was not first created via ``put``
- Use `col_space` argument as minimum column width in DataFrame.to_html (:issue:`2328`)
- Fix tz-aware DatetimeIndex.to_period (:issue:`2232`)
- Fix DataFrame row indexing case with MultiIndex (:issue:`2314`)
- Fix to_excel exporting issues with Timestamp objects in index (:issue:`2294`)
- Fixes assigning scalars and array to hierarchical column chunk (:issue:`1803`)
- Fixed a UnicodeDecodeError with series tidy_repr (:issue:`2225`)
- Fixed issued with duplicate keys in an index (:issue:`2347`, :issue:`2380`)
- Fixed issues re: Hash randomization, default on starting w/ py3.3 (:issue:`2331`)
- Fixed issue with missing attributes after loading a pickled dataframe (:issue:`2431`)
- Fix Timestamp formatting with tzoffset time zone in dateutil 2.1 (:issue:`2443`)
- Fix GroupBy.apply issue when using BinGrouper to do ts binning (:issue:`2300`)
- Fix issues resulting from datetime.datetime columns being converted to
  datetime64 when calling DataFrame.apply. (:issue:`2374`)
- Raise exception when calling to_panel on non uniquely-indexed frame (:issue:`2441`)
- Improved detection of console encoding on IPython zmq frontends (:issue:`2458`)
- Preserve time zone when .appending two time series (:issue:`2260`)
- Box timestamps when calling reset_index on time-zone-aware index rather
  than creating a tz-less datetime64 column (:issue:`2262`)
- Enable searching non-string columns in DataFrame.filter(like=...) (:issue:`2467`)
- Fixed issue with losing nanosecond precision upon conversion to DatetimeIndex(:issue:`2252`)
- Handle timezones in Datetime.normalize (:issue:`2338`)
- Fix test case where dtype specification with endianness causes
  failures on big endian machines (:issue:`2318`)
- Fix plotting bug where upsampling causes data to appear shifted in time (:issue:`2448`)
- Fix ``read_csv`` failure for UTF-16 with BOM and skiprows(:issue:`2298`)
- read_csv with names arg not implicitly setting header=None(:issue:`2459`)
- Unrecognized compression mode causes segfault in read_csv(:issue:`2474`)
- In read_csv, header=0 and passed names should discard first row(:issue:`2269`)
- Correctly route to stdout/stderr in read_table (:issue:`2071`)
- Fix exception when Timestamp.to_datetime is called on a Timestamp with tzoffset (:issue:`2471`)
- Fixed unintentional conversion of datetime64 to long in groupby.first() (:issue:`2133`)
- Union of empty DataFrames now return empty with concatenated index (:issue:`2307`)
- DataFrame.sort_index raises more helpful exception if sorting by column
  with duplicates (:issue:`2488`)
- DataFrame.to_string formatters can be list, too (:issue:`2520`)
- DataFrame.combine_first will always result in the union of the index and
  columns, even if one DataFrame is length-zero (:issue:`2525`)
- Fix several DataFrame.icol/irow with duplicate indices issues (:issue:`2228`, :issue:`2259`)
- Use Series names for column names when using concat with axis=1 (:issue:`2489`)
- Raise Exception if start, end, periods all passed to date_range (:issue:`2538`)
- Fix Panel resampling issue (:issue:`2537`)

pandas 0.9.1
------------

**Release date:** 2012-11-14

New Features
~~~~~~~~~~~~

- Can specify multiple sort orders in DataFrame/Series.sort/sort_index (:issue:`928`)
- New `top` and `bottom` options for handling NAs in rank (:issue:`1508`, :issue:`2159`)
- Add `where` and `mask` functions to DataFrame (:issue:`2109`, :issue:`2151`)
- Add `at_time` and `between_time` functions to DataFrame (:issue:`2149`)
- Add flexible `pow` and `rpow` methods to DataFrame (:issue:`2190`)

API Changes
~~~~~~~~~~~

- Upsampling period index "spans" intervals. Example: annual periods
  upsampled to monthly will span all months in each year
- Period.end_time will yield timestamp at last nanosecond in the interval
  (:issue:`2124`, :issue:`2125`, :issue:`1764`)
- File parsers no longer coerce to float or bool for columns that have custom
  converters specified (:issue:`2184`)

Improvements to existing features
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- Time rule inference for week-of-month (e.g. WOM-2FRI) rules (:issue:`2140`)
- Improve performance of datetime + business day offset with large number of
  offset periods
- Improve HTML display of DataFrame objects with hierarchical columns
- Enable referencing of Excel columns by their column names (:issue:`1936`)
- DataFrame.dot can accept ndarrays (:issue:`2042`)
- Support negative periods in Panel.shift (:issue:`2164`)
- Make .drop(...) work with non-unique indexes (:issue:`2101`)
- Improve performance of Series/DataFrame.diff (re: :issue:`2087`)
- Support unary ~ (__invert__) in DataFrame (:issue:`2110`)
- Turn off pandas-style tick locators and formatters (:issue:`2205`)
- DataFrame[DataFrame] uses DataFrame.where to compute masked frame (:issue:`2230`)

Bug Fixes
~~~~~~~~~

- Fix some duplicate-column DataFrame constructor issues (:issue:`2079`)
- Fix bar plot color cycle issues (:issue:`2082`)
- Fix off-center grid for stacked bar plots (:issue:`2157`)
- Fix plotting bug if inferred frequency is offset with N > 1 (:issue:`2126`)
- Implement comparisons on date offsets with fixed delta (:issue:`2078`)
- Handle inf/-inf correctly in read_* parser functions (:issue:`2041`)
- Fix matplotlib unicode interaction bug
- Make WLS r-squared match statsmodels 0.5.0 fixed value
- Fix zero-trimming DataFrame formatting bug
- Correctly compute/box datetime64 min/max values from Series.min/max (:issue:`2083`)
- Fix unstacking edge case with unrepresented groups (:issue:`2100`)
- Fix Series.str failures when using pipe pattern '|' (:issue:`2119`)
- Fix pretty-printing of dict entries in Series, DataFrame (:issue:`2144`)
- Cast other datetime64 values to nanoseconds in DataFrame ctor (:issue:`2095`)
- Alias Timestamp.astimezone to tz_convert, so will yield Timestamp (:issue:`2060`)
- Fix timedelta64 formatting from Series (:issue:`2165`, :issue:`2146`)
- Handle None values gracefully in dict passed to Panel constructor (:issue:`2075`)
- Box datetime64 values as Timestamp objects in Series/DataFrame.iget (:issue:`2148`)
- Fix Timestamp indexing bug in DatetimeIndex.insert (:issue:`2155`)
- Use index name(s) (if any) in DataFrame.to_records (:issue:`2161`)
- Don't lose index names in Panel.to_frame/DataFrame.to_panel (:issue:`2163`)
- Work around length-0 boolean indexing NumPy bug (:issue:`2096`)
- Fix partial integer indexing bug in DataFrame.xs (:issue:`2107`)
- Fix variety of cut/qcut string-bin formatting bugs (:issue:`1978`, :issue:`1979`)
- Raise Exception when xs view not possible of MultiIndex'd DataFrame (:issue:`2117`)
- Fix groupby(...).first() issue with datetime64 (:issue:`2133`)
- Better floating point error robustness in some rolling_* functions
  (:issue:`2114`, :issue:`2527`)
- Fix ewma NA handling in the middle of Series (:issue:`2128`)
- Fix numerical precision issues in diff with integer data (:issue:`2087`)
- Fix bug in MultiIndex.__getitem__ with NA values (:issue:`2008`)
- Fix DataFrame.from_records dict-arg bug when passing columns (:issue:`2179`)
- Fix Series and DataFrame.diff for integer dtypes (:issue:`2087`, :issue:`2174`)
- Fix bug when taking intersection of DatetimeIndex with empty index (:issue:`2129`)
- Pass through timezone information when calling DataFrame.align (:issue:`2127`)
- Properly sort when joining on datetime64 values (:issue:`2196`)
- Fix indexing bug in which False/True were being coerced to 0/1 (:issue:`2199`)
- Many unicode formatting fixes (:issue:`2201`)
- Fix improper MultiIndex conversion issue when assigning
  e.g. DataFrame.index (:issue:`2200`)
- Fix conversion of mixed-type DataFrame to ndarray with dup columns (:issue:`2236`)
- Fix duplicate columns issue (:issue:`2218`, :issue:`2219`)
- Fix SparseSeries.__pow__ issue with NA input (:issue:`2220`)
- Fix icol with integer sequence failure (:issue:`2228`)
- Fixed resampling tz-aware time series issue (:issue:`2245`)
- SparseDataFrame.icol was not returning SparseSeries (:issue:`2227`, :issue:`2229`)
- Enable ExcelWriter to handle PeriodIndex (:issue:`2240`)
- Fix issue constructing DataFrame from empty Series with name (:issue:`2234`)
- Use console-width detection in interactive sessions only (:issue:`1610`)
- Fix parallel_coordinates legend bug with mpl 1.2.0 (:issue:`2237`)
- Make tz_localize work in corner case of empty Series (:issue:`2248`)

pandas 0.9.0
------------

**Release date:** 10/7/2012

New Features
~~~~~~~~~~~~

- Add ``str.encode`` and ``str.decode`` to Series (:issue:`1706`)
- Add `to_latex` method to DataFrame (:issue:`1735`)
- Add convenient expanding window equivalents of all rolling_* ops (:issue:`1785`)
- Add Options class to pandas.io.data for fetching options data from Yahoo!
  Finance (:issue:`1748`, :issue:`1739`)
- Recognize and convert more boolean values in file parsing (Yes, No, TRUE,
  FALSE, variants thereof) (:issue:`1691`, :issue:`1295`)
- Add Panel.update method, analogous to DataFrame.update (:issue:`1999`, :issue:`1988`)

Improvements to existing features
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- Proper handling of NA values in merge operations (:issue:`1990`)
- Add ``flags`` option for ``re.compile`` in some Series.str methods (:issue:`1659`)
- Parsing of UTC date strings in read_* functions (:issue:`1693`)
- Handle generator input to Series (:issue:`1679`)
- Add `na_action='ignore'` to Series.map to quietly propagate NAs (:issue:`1661`)
- Add args/kwds options to Series.apply (:issue:`1829`)
- Add inplace option to Series/DataFrame.reset_index (:issue:`1797`)
- Add ``level`` parameter to ``Series.reset_index``
- Add quoting option for DataFrame.to_csv (:issue:`1902`)
- Indicate long column value truncation in DataFrame output with ... (:issue:`1854`)
- DataFrame.dot will not do data alignment, and also work with Series (:issue:`1915`)
- Add ``na`` option for missing data handling in some vectorized string
  methods (:issue:`1689`)
- If index_label=False in DataFrame.to_csv, do not print fields/commas in the
  text output. Results in easier importing into R (:issue:`1583`)
- Can pass tuple/list of axes to DataFrame.dropna to simplify repeated calls
  (dropping both columns and rows) (:issue:`924`)
- Improve DataFrame.to_html output for hierarchically-indexed rows (do not
  repeat levels) (:issue:`1929`)
- TimeSeries.between_time can now select times across midnight (:issue:`1871`)
- Enable `skip_footer` parameter in `ExcelFile.parse` (:issue:`1843`)

API Changes
~~~~~~~~~~~

- Change default header names in read_* functions to more Pythonic X0, X1,
  etc. instead of X.1, X.2. (:issue:`2000`)
- Deprecated ``day_of_year`` API removed from PeriodIndex, use ``dayofyear``
  (:issue:`1723`)
- Don't modify NumPy suppress printoption at import time
- The internal HDF5 data arrangement for DataFrames has been
  transposed. Legacy files will still be readable by HDFStore (:issue:`1834`, :issue:`1824`)
- Legacy cruft removed: pandas.stats.misc.quantileTS
- Use ISO8601 format for Period repr: monthly, daily, and on down (:issue:`1776`)
- Empty DataFrame columns are now created as object dtype. This will prevent
  a class of TypeErrors that was occurring in code where the dtype of a
  column would depend on the presence of data or not (e.g. a SQL query having
  results) (:issue:`1783`)
- Setting parts of DataFrame/Panel using ix now aligns input Series/DataFrame
  (:issue:`1630`)
- `first` and `last` methods in `GroupBy` no longer drop non-numeric columns
  (:issue:`1809`)
- Resolved inconsistencies in specifying custom NA values in text parser.
  `na_values` of type dict no longer override default NAs unless
  `keep_default_na` is set to false explicitly (:issue:`1657`)
- Enable `skipfooter` parameter in text parsers as an alias for `skip_footer`

Bug Fixes
~~~~~~~~~

- Perform arithmetic column-by-column in mixed-type DataFrame to avoid type
  upcasting issues. Caused downstream DataFrame.diff bug (:issue:`1896`)
- Fix matplotlib auto-color assignment when no custom spectrum passed. Also
  respect passed color keyword argument (:issue:`1711`)
- Fix resampling logical error with closed='left' (:issue:`1726`)
- Fix critical DatetimeIndex.union bugs (:issue:`1730`, :issue:`1719`, :issue:`1745`, :issue:`1702`, :issue:`1753`)
- Fix critical DatetimeIndex.intersection bug with unanchored offsets (:issue:`1708`)
- Fix MM-YYYY time series indexing case (:issue:`1672`)
- Fix case where Categorical group key was not being passed into index in
  GroupBy result (:issue:`1701`)
- Handle Ellipsis in Series.__getitem__/__setitem__ (:issue:`1721`)
- Fix some bugs with handling datetime64 scalars of other units in NumPy 1.6
  and 1.7 (:issue:`1717`)
- Fix performance issue in MultiIndex.format (:issue:`1746`)
- Fixed GroupBy bugs interacting with DatetimeIndex asof / map methods (:issue:`1677`)
- Handle factors with NAs in pandas.rpy (:issue:`1615`)
- Fix statsmodels import in pandas.stats.var (:issue:`1734`)
- Fix DataFrame repr/info summary with non-unique columns (:issue:`1700`)
- Fix Series.iget_value for non-unique indexes (:issue:`1694`)
- Don't lose tzinfo when passing DatetimeIndex as DataFrame column (:issue:`1682`)
- Fix tz conversion with time zones that haven't had any DST transitions since
  first date in the array (:issue:`1673`)
- Fix field access with UTC->local conversion on unsorted arrays (:issue:`1756`)
- Fix isnull handling of array-like (list) inputs (:issue:`1755`)
- Fix regression in handling of Series in Series constructor (:issue:`1671`)
- Fix comparison of Int64Index with DatetimeIndex (:issue:`1681`)
- Fix min_periods handling in new rolling_max/min at array start (:issue:`1695`)
- Fix errors with how='median' and generic NumPy resampling in some cases
  caused by SeriesBinGrouper (:issue:`1648`, :issue:`1688`)
- When grouping by level, exclude unobserved levels (:issue:`1697`)
- Don't lose tzinfo in DatetimeIndex when shifting by different offset (:issue:`1683`)
- Hack to support storing data with a zero-length axis in HDFStore (:issue:`1707`)
- Fix DatetimeIndex tz-aware range generation issue (:issue:`1674`)
- Fix method='time' interpolation with intraday data (:issue:`1698`)
- Don't plot all-NA DataFrame columns as zeros (:issue:`1696`)
- Fix bug in scatter_plot with by option (:issue:`1716`)
- Fix performance problem in infer_freq with lots of non-unique stamps (:issue:`1686`)
- Fix handling of PeriodIndex as argument to create MultiIndex (:issue:`1705`)
- Fix re: unicode MultiIndex level names in Series/DataFrame repr (:issue:`1736`)
- Handle PeriodIndex in to_datetime instance method (:issue:`1703`)
- Support StaticTzInfo in DatetimeIndex infrastructure (:issue:`1692`)
- Allow MultiIndex setops with length-0 other type indexes (:issue:`1727`)
- Fix handling of DatetimeIndex in DataFrame.to_records (:issue:`1720`)
- Fix handling of general objects in isnull on which bool(...) fails (:issue:`1749`)
- Fix .ix indexing with MultiIndex ambiguity (:issue:`1678`)
- Fix .ix setting logic error with non-unique MultiIndex (:issue:`1750`)
- Basic indexing now works on MultiIndex with > 1000000 elements, regression
  from earlier version of pandas (:issue:`1757`)
- Handle non-float64 dtypes in fast DataFrame.corr/cov code paths (:issue:`1761`)
- Fix DatetimeIndex.isin to function properly (:issue:`1763`)
- Fix conversion of array of tz-aware datetime.datetime to DatetimeIndex with
  right time zone (:issue:`1777`)
- Fix DST issues with generating anchored date ranges (:issue:`1778`)
- Fix issue calling sort on result of Series.unique (:issue:`1807`)
- Fix numerical issue leading to square root of negative number in
  rolling_std (:issue:`1840`)
- Let Series.str.split accept no arguments (like str.split) (:issue:`1859`)
- Allow user to have dateutil 2.1 installed on a Python 2 system (:issue:`1851`)
- Catch ImportError less aggressively in pandas/__init__.py (:issue:`1845`)
- Fix pip source installation bug when installing from GitHub (:issue:`1805`)
- Fix error when window size > array size in rolling_apply (:issue:`1850`)
- Fix pip source installation issues via SSH from GitHub
- Fix OLS.summary when column is a tuple (:issue:`1837`)
- Fix bug in __doc__ patching when -OO passed to interpreter
  (:issue:`1792` :issue:`1741` :issue:`1774`)
- Fix unicode console encoding issue in IPython notebook (:issue:`1782`, :issue:`1768`)
- Fix unicode formatting issue with Series.name (:issue:`1782`)
- Fix bug in DataFrame.duplicated with datetime64 columns (:issue:`1833`)
- Fix bug in Panel internals resulting in error when doing fillna after
  truncate not changing size of panel (:issue:`1823`)
- Prevent segfault due to MultiIndex not being supported in HDFStore table
  format (:issue:`1848`)
- Fix UnboundLocalError in Panel.__setitem__ and add better error (:issue:`1826`)
- Fix to_csv issues with list of string entries. Isnull works on list of
  strings now too (:issue:`1791`)
- Fix Timestamp comparisons with datetime values outside the nanosecond range
  (1677-2262)
- Revert to prior behavior of normalize_date with datetime.date objects
  (return datetime)
- Fix broken interaction between np.nansum and Series.any/all
- Fix bug with multiple column date parsers (:issue:`1866`)
- DatetimeIndex.union(Int64Index) was broken
- Make plot x vs y interface consistent with integer indexing (:issue:`1842`)
- set_index inplace modified data even if unique check fails (:issue:`1831`)
- Only use Q-OCT/NOV/DEC in quarterly frequency inference (:issue:`1789`)
- Upcast to dtype=object when unstacking boolean DataFrame (:issue:`1820`)
- Fix float64/float32 merging bug (:issue:`1849`)
- Fixes to Period.start_time for non-daily frequencies (:issue:`1857`)
- Fix failure when converter used on index_col in read_csv (:issue:`1835`)
- Implement PeriodIndex.append so that pandas.concat works correctly (:issue:`1815`)
- Avoid Cython out-of-bounds access causing segfault sometimes in pad_2d,
  backfill_2d
- Fix resampling error with intraday times and anchored target time (like
  AS-DEC) (:issue:`1772`)
- Fix .ix indexing bugs with mixed-integer indexes (:issue:`1799`)
- Respect passed color keyword argument in Series.plot (:issue:`1890`)
- Fix rolling_min/max when the window is larger than the size of the input
  array. Check other malformed inputs (:issue:`1899`, :issue:`1897`)
- Rolling variance / standard deviation with only a single observation in
  window (:issue:`1884`)
- Fix unicode sheet name failure in to_excel (:issue:`1828`)
- Override DatetimeIndex.min/max to return Timestamp objects (:issue:`1895`)
- Fix column name formatting issue in length-truncated column (:issue:`1906`)
- Fix broken handling of copying Index metadata to new instances created by
  view(...) calls inside the NumPy infrastructure
- Support datetime.date again in DateOffset.rollback/rollforward
- Raise Exception if set passed to Series constructor (:issue:`1913`)
- Add TypeError when appending HDFStore table w/ wrong index type (:issue:`1881`)
- Don't raise exception on empty inputs in EW functions (e.g. ewma) (:issue:`1900`)
- Make asof work correctly with PeriodIndex (:issue:`1883`)
- Fix extlinks in doc build
- Fill boolean DataFrame with NaN when calling shift (:issue:`1814`)
- Fix setuptools bug causing pip not to Cythonize .pyx files sometimes
- Fix negative integer indexing regression in .ix from 0.7.x (:issue:`1888`)
- Fix error while retrieving timezone and utc offset from subclasses of
  datetime.tzinfo without .zone and ._utcoffset attributes (:issue:`1922`)
- Fix DataFrame formatting of small, non-zero FP numbers (:issue:`1911`)
- Various fixes by upcasting of date -> datetime (:issue:`1395`)
- Raise better exception when passing multiple functions with the same name,
  such as lambdas, to GroupBy.aggregate
- Fix DataFrame.apply with axis=1 on a non-unique index (:issue:`1878`)
- Proper handling of Index subclasses in pandas.unique (:issue:`1759`)
- Set index names in DataFrame.from_records (:issue:`1744`)
- Fix time series indexing error with duplicates, under and over hash table
  size cutoff (:issue:`1821`)
- Handle list keys in addition to tuples in DataFrame.xs when
  partial-indexing a hierarchically-indexed DataFrame (:issue:`1796`)
- Support multiple column selection in DataFrame.__getitem__ with duplicate
  columns (:issue:`1943`)
- Fix time zone localization bug causing improper fields (e.g. hours) in time
  zones that have not had a UTC transition in a long time (:issue:`1946`)
- Fix errors when parsing and working with fixed offset timezones
  (:issue:`1922`, :issue:`1928`)
- Fix text parser bug when handling UTC datetime objects generated by
  dateutil (:issue:`1693`)
- Fix plotting bug when 'B' is the inferred frequency but index actually
  contains weekends (:issue:`1668`, :issue:`1669`)
- Fix plot styling bugs (:issue:`1666`, :issue:`1665`, :issue:`1658`)
- Fix plotting bug with index/columns with unicode (:issue:`1685`)
- Fix DataFrame constructor bug when passed Series with datetime64 dtype
  in a dict (:issue:`1680`)
- Fixed regression in generating DatetimeIndex using timezone aware
  datetime.datetime (:issue:`1676`)
- Fix DataFrame bug when printing concatenated DataFrames with duplicated
  columns (:issue:`1675`)
- Fixed bug when plotting time series with multiple intraday frequencies
  (:issue:`1732`)
- Fix bug in DataFrame.duplicated to enable iterables other than list-types
  as input argument (:issue:`1773`)
- Fix resample bug when passed list of lambdas as `how` argument (:issue:`1808`)
- Repr fix for MultiIndex level with all NAs (:issue:`1971`)
- Fix PeriodIndex slicing bug when slice start/end are out-of-bounds (:issue:`1977`)
- Fix read_table bug when parsing unicode (:issue:`1975`)
- Fix BlockManager.iget bug when dealing with non-unique MultiIndex as columns
  (:issue:`1970`)
- Fix reset_index bug if both drop and level are specified (:issue:`1957`)
- Work around unsafe NumPy object->int casting with Cython function (:issue:`1987`)
- Fix datetime64 formatting bug in DataFrame.to_csv (:issue:`1993`)
- Default start date in pandas.io.data to 1/1/2000 as the docs say (:issue:`2011`)

pandas 0.8.1
------------

**Release date:** July 22, 2012

New Features
~~~~~~~~~~~~

- Add vectorized, NA-friendly string methods to Series (:issue:`1621`, :issue:`620`)
- Can pass dict of per-column line styles to DataFrame.plot (:issue:`1559`)
- Selective plotting to secondary y-axis on same subplot (:issue:`1640`)
- Add new ``bootstrap_plot`` plot function
- Add new ``parallel_coordinates`` plot function (:issue:`1488`)
- Add ``radviz`` plot function (:issue:`1566`)
- Add ``multi_sparse`` option to ``set_printoptions`` to modify display of
  hierarchical indexes (:issue:`1538`)
- Add ``dropna`` method to Panel (:issue:`171`)

Improvements to existing features
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- Use moving min/max algorithms from Bottleneck in rolling_min/rolling_max
  for > 100x speedup. (:issue:`1504`, :issue:`50`)
- Add Cython group median method for >15x speedup (:issue:`1358`)
- Drastically improve ``to_datetime`` performance on ISO8601 datetime strings
  (with no time zones) (:issue:`1571`)
- Improve single-key groupby performance on large data sets, accelerate use of
  groupby with a Categorical variable
- Add ability to append hierarchical index levels with ``set_index`` and to
  drop single levels with ``reset_index`` (:issue:`1569`, :issue:`1577`)
- Always apply passed functions in ``resample``, even if upsampling (:issue:`1596`)
- Avoid unnecessary copies in DataFrame constructor with explicit dtype (:issue:`1572`)
- Cleaner DatetimeIndex string representation with 1 or 2 elements (:issue:`1611`)
- Improve performance of array-of-Period to PeriodIndex, convert such arrays
  to PeriodIndex inside Index (:issue:`1215`)
- More informative string representation for weekly Period objects (:issue:`1503`)
- Accelerate 3-axis multi data selection from homogeneous Panel (:issue:`979`)
- Add ``adjust`` option to ewma to disable adjustment factor (:issue:`1584`)
- Add new matplotlib converters for high frequency time series plotting (:issue:`1599`)
- Handling of tz-aware datetime.datetime objects in to_datetime; raise
  Exception unless utc=True given (:issue:`1581`)

Bug Fixes
~~~~~~~~~

- Fix NA handling in DataFrame.to_panel (:issue:`1582`)
- Handle TypeError issues inside PyObject_RichCompareBool calls in khash
  (:issue:`1318`)
- Fix resampling bug to lower case daily frequency (:issue:`1588`)
- Fix kendall/spearman DataFrame.corr bug with no overlap (:issue:`1595`)
- Fix bug in DataFrame.set_index (:issue:`1592`)
- Don't ignore axes in boxplot if by specified (:issue:`1565`)
- Fix Panel .ix indexing with integers bug (:issue:`1603`)
- Fix Partial indexing bugs (years, months, ...) with PeriodIndex (:issue:`1601`)
- Fix MultiIndex console formatting issue (:issue:`1606`)
- Unordered index with duplicates doesn't yield scalar location for single
  entry (:issue:`1586`)
- Fix resampling of tz-aware time series with "anchored" freq (:issue:`1591`)
- Fix DataFrame.rank error on integer data (:issue:`1589`)
- Selection of multiple SparseDataFrame columns by list in __getitem__ (:issue:`1585`)
- Override Index.tolist for compatibility with MultiIndex (:issue:`1576`)
- Fix hierarchical summing bug with MultiIndex of length 1 (:issue:`1568`)
- Work around numpy.concatenate use/bug in Series.set_value (:issue:`1561`)
- Ensure Series/DataFrame are sorted before resampling (:issue:`1580`)
- Fix unhandled IndexError when indexing very large time series (:issue:`1562`)
- Fix DatetimeIndex intersection logic error with irregular indexes (:issue:`1551`)
- Fix unit test errors on Python 3 (:issue:`1550`)
- Fix .ix indexing bugs in duplicate DataFrame index (:issue:`1201`)
- Better handle errors with non-existing objects in HDFStore (:issue:`1254`)
- Don't copy int64 array data in DatetimeIndex when copy=False (:issue:`1624`)
- Fix resampling of conforming periods quarterly to annual (:issue:`1622`)
- Don't lose index name on resampling (:issue:`1631`)
- Support python-dateutil version 2.1 (:issue:`1637`)
- Fix broken scatter_matrix axis labeling, esp. with time series (:issue:`1625`)
- Fix cases where extra keywords weren't being passed on to matplotlib from
  Series.plot (:issue:`1636`)
- Fix BusinessMonthBegin logic for dates before 1st bday of month (:issue:`1645`)
- Ensure string alias converted (valid in DatetimeIndex.get_loc) in
  DataFrame.xs / __getitem__ (:issue:`1644`)
- Fix use of string alias timestamps with tz-aware time series (:issue:`1647`)
- Fix Series.max/min and Series.describe on len-0 series (:issue:`1650`)
- Handle None values in dict passed to concat (:issue:`1649`)
- Fix Series.interpolate with method='values' and DatetimeIndex (:issue:`1646`)
- Fix IndexError in left merges on a DataFrame with 0-length (:issue:`1628`)
- Fix DataFrame column width display with UTF-8 encoded characters (:issue:`1620`)
- Handle case in pandas.io.data.get_data_yahoo where Yahoo! returns duplicate
  dates for most recent business day
- Avoid downsampling when plotting mixed frequencies on the same subplot (:issue:`1619`)
- Fix read_csv bug when reading a single line (:issue:`1553`)
- Fix bug in C code causing monthly periods prior to December 1969 to be off (:issue:`1570`)

pandas 0.8.0
------------

**Release date:** 6/29/2012

New Features
~~~~~~~~~~~~

- New unified DatetimeIndex class for nanosecond-level timestamp data
- New Timestamp datetime.datetime subclass with easy time zone conversions,
  and support for nanoseconds
- New PeriodIndex class for timespans, calendar logic, and Period scalar object
- High performance resampling of timestamp and period data. New `resample`
  method of all pandas data structures
- New frequency names plus shortcut string aliases like '15h', '1h30min'
- Time series string indexing shorthand (:issue:`222`)
- Add week, dayofyear array and other timestamp array-valued field accessor
  functions to DatetimeIndex
- Add GroupBy.prod optimized aggregation function and 'prod' fast time series
  conversion method (:issue:`1018`)
- Implement robust frequency inference function and `inferred_freq` attribute
  on DatetimeIndex (:issue:`391`)
- New ``tz_convert`` and ``tz_localize`` methods in Series / DataFrame
- Convert DatetimeIndexes to UTC if time zones are different in join/setops
  (:issue:`864`)
- Add limit argument for forward/backward filling to reindex, fillna,
  etc. (:issue:`825` and others)
- Add support for indexes (dates or otherwise) with duplicates and common
  sense indexing/selection functionality
- Series/DataFrame.update methods, in-place variant of combine_first (:issue:`961`)
- Add ``match`` function to API (:issue:`502`)
- Add Cython-optimized first, last, min, max, prod functions to GroupBy (:issue:`994`,
  :issue:`1043`)
- Dates can be split across multiple columns (:issue:`1227`, :issue:`1186`)
- Add experimental support for converting pandas DataFrame to R data.frame
  via rpy2 (:issue:`350`, :issue:`1212`)
- Can pass list of (name, function) to GroupBy.aggregate to get aggregates in
  a particular order (:issue:`610`)
- Can pass dicts with lists of functions or dicts to GroupBy aggregate to do
  much more flexible multiple function aggregation (:issue:`642`, :issue:`610`)
- New ordered_merge functions for merging DataFrames with ordered
  data. Also supports group-wise merging for panel data (:issue:`813`)
- Add keys() method to DataFrame
- Add flexible replace method for replacing potentially values to Series and
  DataFrame (:issue:`929`, :issue:`1241`)
- Add 'kde' plot kind for Series/DataFrame.plot (:issue:`1059`)
- More flexible multiple function aggregation with GroupBy
- Add pct_change function to Series/DataFrame
- Add option to interpolate by Index values in Series.interpolate (:issue:`1206`)
- Add ``max_colwidth`` option for DataFrame, defaulting to 50
- Conversion of DataFrame through rpy2 to R data.frame (:issue:`1282`, )
- Add keys() method on DataFrame (:issue:`1240`)
- Add new ``match`` function to API (similar to R) (:issue:`502`)
- Add dayfirst option to parsers (:issue:`854`)
- Add ``method`` argument to ``align`` method for forward/backward filling
  (:issue:`216`)
- Add Panel.transpose method for rearranging axes (:issue:`695`)
- Add new ``cut`` function (patterned after R) for discretizing data into
  equal range-length bins or arbitrary breaks of your choosing (:issue:`415`)
- Add new ``qcut`` for cutting with quantiles (:issue:`1378`)
- Add ``value_counts`` top level array method (:issue:`1392`)
- Added Andrews curves plot type (:issue:`1325`)
- Add lag plot (:issue:`1440`)
- Add autocorrelation_plot (:issue:`1425`)
- Add support for tox and Travis CI (:issue:`1382`)
- Add support for Categorical use in GroupBy (:issue:`292`)
- Add ``any`` and ``all`` methods to DataFrame (:issue:`1416`)
- Add ``secondary_y`` option to Series.plot
- Add experimental ``lreshape`` function for reshaping wide to long

Improvements to existing features
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- Switch to klib/khash-based hash tables in Index classes for better
  performance in many cases and lower memory footprint
- Shipping some functions from scipy.stats to reduce dependency,
  e.g. Series.describe and DataFrame.describe (:issue:`1092`)
- Can create MultiIndex by passing list of lists or list of arrays to Series,
  DataFrame constructor, etc. (:issue:`831`)
- Can pass arrays in addition to column names to DataFrame.set_index (:issue:`402`)
- Improve the speed of "square" reindexing of homogeneous DataFrame objects
  by significant margin (:issue:`836`)
- Handle more dtypes when passed MaskedArrays in DataFrame constructor (:issue:`406`)
- Improved performance of join operations on integer keys (:issue:`682`)
- Can pass multiple columns to GroupBy object, e.g. grouped[[col1, col2]] to
  only aggregate a subset of the value columns (:issue:`383`)
- Add histogram / kde plot options for scatter_matrix diagonals (:issue:`1237`)
- Add inplace option to Series/DataFrame.rename and sort_index,
  DataFrame.drop_duplicates (:issue:`805`, :issue:`207`)
- More helpful error message when nothing passed to Series.reindex (:issue:`1267`)
- Can mix array and scalars as dict-value inputs to DataFrame ctor (:issue:`1329`)
- Use DataFrame columns' name for legend title in plots
- Preserve frequency in DatetimeIndex when possible in boolean indexing
  operations
- Promote datetime.date values in data alignment operations (:issue:`867`)
- Add ``order`` method to Index classes (:issue:`1028`)
- Avoid hash table creation in large monotonic hash table indexes (:issue:`1160`)
- Store time zones in HDFStore (:issue:`1232`)
- Enable storage of sparse data structures in HDFStore (:issue:`85`)
- Enable Series.asof to work with arrays of timestamp inputs
- Cython implementation of DataFrame.corr speeds up by > 100x (:issue:`1349`, :issue:`1354`)
- Exclude "nuisance" columns automatically in GroupBy.transform (:issue:`1364`)
- Support functions-as-strings in GroupBy.transform (:issue:`1362`)
- Use index name as xlabel/ylabel in plots (:issue:`1415`)
- Add ``convert_dtype`` option to Series.apply to be able to leave data as
  dtype=object (:issue:`1414`)
- Can specify all index level names in concat (:issue:`1419`)
- Add ``dialect`` keyword to parsers for quoting conventions (:issue:`1363`)
- Enable DataFrame[bool_DataFrame] += value (:issue:`1366`)
- Add ``retries`` argument to ``get_data_yahoo`` to try to prevent Yahoo! API
  404s (:issue:`826`)
- Improve performance of reshaping by using O(N) categorical sorting
- Series names will be used for index of DataFrame if no index passed (:issue:`1494`)
- Header argument in DataFrame.to_csv can accept a list of column names to
  use instead of the object's columns (:issue:`921`)
- Add ``raise_conflict`` argument to DataFrame.update (:issue:`1526`)
- Support file-like objects in ExcelFile (:issue:`1529`)

API Changes
~~~~~~~~~~~

- Rename `pandas._tseries` to `pandas.lib`
- Rename Factor to Categorical and add improvements. Numerous Categorical bug
  fixes
- Frequency name overhaul, WEEKDAY/EOM and rules with @
  deprecated. get_legacy_offset_name backwards compatibility function added
- Raise ValueError in DataFrame.__nonzero__, so "if df" no longer works
  (:issue:`1073`)
- Change BDay (business day) to not normalize dates by default (:issue:`506`)
- Remove deprecated DataMatrix name
- Default merge suffixes for overlap now have underscores instead of periods
  to facilitate tab completion, etc. (:issue:`1239`)
- Deprecation of offset, time_rule timeRule parameters throughout code base
- Series.append and DataFrame.append no longer check for duplicate indexes
  by default, add verify_integrity parameter (:issue:`1394`)
- Refactor Factor class, old constructor moved to Factor.from_array
- Modified internals of MultiIndex to use less memory (no longer represented
  as array of tuples) internally, speed up construction time and many methods
  which construct intermediate hierarchical indexes (:issue:`1467`)

Bug Fixes
~~~~~~~~~

- Fix OverflowError from storing pre-1970 dates in HDFStore by switching to
  datetime64 (:issue:`179`)
- Fix logical error with February leap year end in YearEnd offset
- Series([False, nan]) was getting casted to float64 (:issue:`1074`)
- Fix binary operations between boolean Series and object Series with
  booleans and NAs (:issue:`1074`, :issue:`1079`)
- Couldn't assign whole array to column in mixed-type DataFrame via .ix
  (:issue:`1142`)
- Fix label slicing issues with float index values (:issue:`1167`)
- Fix segfault caused by empty groups passed to groupby (:issue:`1048`)
- Fix occasionally misbehaved reindexing in the presence of NaN labels (:issue:`522`)
- Fix imprecise logic causing weird Series results from .apply (:issue:`1183`)
- Unstack multiple levels in one shot, avoiding empty columns in some
  cases. Fix pivot table bug (:issue:`1181`)
- Fix formatting of MultiIndex on Series/DataFrame when index name coincides
  with label (:issue:`1217`)
- Handle Excel 2003 #N/A as NaN from xlrd (:issue:`1213`, :issue:`1225`)
- Fix timestamp locale-related deserialization issues with HDFStore by moving
  to datetime64 representation (:issue:`1081`, :issue:`809`)
- Fix DataFrame.duplicated/drop_duplicates NA value handling (:issue:`557`)
- Actually raise exceptions in fast reducer (:issue:`1243`)
- Fix various timezone-handling bugs from 0.7.3 (:issue:`969`)
- GroupBy on level=0 discarded index name (:issue:`1313`)
- Better error message with unmergeable DataFrames (:issue:`1307`)
- Series.__repr__ alignment fix with unicode index values (:issue:`1279`)
- Better error message if nothing passed to reindex (:issue:`1267`)
- More robust NA handling in DataFrame.drop_duplicates (:issue:`557`)
- Resolve locale-based and pre-epoch HDF5 timestamp deserialization issues
  (:issue:`973`, :issue:`1081`, :issue:`179`)
- Implement Series.repeat (:issue:`1229`)
- Fix indexing with namedtuple and other tuple subclasses (:issue:`1026`)
- Fix float64 slicing bug (:issue:`1167`)
- Parsing integers with commas (:issue:`796`)
- Fix groupby improper data type when group consists of one value (:issue:`1065`)
- Fix negative variance possibility in nanvar resulting from floating point
  error (:issue:`1090`)
- Consistently set name on groupby pieces (:issue:`184`)
- Treat dict return values as Series in GroupBy.apply (:issue:`823`)
- Respect column selection for DataFrame in GroupBy.transform (:issue:`1365`)
- Fix MultiIndex partial indexing bug (:issue:`1352`)
- Enable assignment of rows in mixed-type DataFrame via .ix (:issue:`1432`)
- Reset index mapping when grouping Series in Cython (:issue:`1423`)
- Fix outer/inner DataFrame.join with non-unique indexes (:issue:`1421`)
- Fix MultiIndex groupby bugs with empty lower levels (:issue:`1401`)
- Calling fillna with a Series will have same behavior as with dict (:issue:`1486`)
- SparseSeries reduction bug (:issue:`1375`)
- Fix unicode serialization issue in HDFStore (:issue:`1361`)
- Pass keywords to pyplot.boxplot in DataFrame.boxplot (:issue:`1493`)
- Bug fixes in MonthBegin (:issue:`1483`)
- Preserve MultiIndex names in drop (:issue:`1513`)
- Fix Panel DataFrame slice-assignment bug (:issue:`1533`)
- Don't use locals() in read_* functions (:issue:`1547`)

pandas 0.7.3
------------

**Release date:** April 12, 2012

New Features
~~~~~~~~~~~~

- Support for non-unique indexes: indexing and selection, many-to-one and
  many-to-many joins (:issue:`1306`)
- Added fixed-width file reader, read_fwf (:issue:`952`)
- Add group_keys argument to groupby to not add group names to MultiIndex in
  result of apply (:issue:`938`)
- DataFrame can now accept non-integer label slicing (:issue:`946`). Previously
  only DataFrame.ix was able to do so.
- DataFrame.apply now retains name attributes on Series objects (:issue:`983`)
- Numeric DataFrame comparisons with non-numeric values now raises proper
  TypeError (:issue:`943`). Previously raise "PandasError: DataFrame constructor
  not properly called!"
- Add ``kurt`` methods to Series and DataFrame (:issue:`964`)
- Can pass dict of column -> list/set NA values for text parsers (:issue:`754`)
- Allows users specified NA values in text parsers (:issue:`754`)
- Parsers checks for openpyxl dependency and raises ImportError if not found
  (:issue:`1007`)
- New factory function to create HDFStore objects that can be used in a with
  statement so users do not have to explicitly call HDFStore.close (:issue:`1005`)
- pivot_table is now more flexible with same parameters as groupby (:issue:`941`)
- Added stacked bar plots (:issue:`987`)
- scatter_matrix method in pandas/tools/plotting.py (:issue:`935`)
- DataFrame.boxplot returns plot results for ex-post styling (:issue:`985`)
- Short version number accessible as pandas.version.short_version (:issue:`930`)
- Additional documentation in panel.to_frame (:issue:`942`)
- More informative Series.apply docstring regarding element-wise apply
  (:issue:`977`)
- Notes on rpy2 installation (:issue:`1006`)
- Add rotation and font size options to hist method (:issue:`1012`)
- Use exogenous / X variable index in result of OLS.y_predict. Add
  OLS.predict method (:issue:`1027`, :issue:`1008`)

API Changes
~~~~~~~~~~~

- Calling apply on grouped Series, e.g. describe(), will no longer yield
  DataFrame by default. Will have to call unstack() to get prior behavior
- NA handling in non-numeric comparisons has been tightened up (:issue:`933`, :issue:`953`)
- No longer assign dummy names key_0, key_1, etc. to groupby index (:issue:`1291`)

Bug Fixes
~~~~~~~~~

- Fix logic error when selecting part of a row in a DataFrame with a
  MultiIndex index (:issue:`1013`)
- Series comparison with Series of differing length causes crash (:issue:`1016`).
- Fix bug in indexing when selecting section of hierarchically-indexed row
  (:issue:`1013`)
- DataFrame.plot(logy=True) has no effect (:issue:`1011`).
- Broken arithmetic operations between SparsePanel-Panel (:issue:`1015`)
- Unicode repr issues in MultiIndex with non-ASCII characters (:issue:`1010`)
- DataFrame.lookup() returns inconsistent results if exact match not present
  (:issue:`1001`)
- DataFrame arithmetic operations not treating None as NA (:issue:`992`)
- DataFrameGroupBy.apply returns incorrect result (:issue:`991`)
- Series.reshape returns incorrect result for multiple dimensions (:issue:`989`)
- Series.std and Series.var ignores ddof parameter (:issue:`934`)
- DataFrame.append loses index names (:issue:`980`)
- DataFrame.plot(kind='bar') ignores color argument (:issue:`958`)
- Inconsistent Index comparison results (:issue:`948`)
- Improper int dtype DataFrame construction from data with NaN (:issue:`846`)
- Removes default 'result' name in groupby results (:issue:`995`)
- DataFrame.from_records no longer mutate input columns (:issue:`975`)
- Use Index name when grouping by it (:issue:`1313`)

pandas 0.7.2
------------

**Release date:** March 16, 2012

New Features
~~~~~~~~~~~~

- Add additional tie-breaking methods in DataFrame.rank (:issue:`874`)
- Add ascending parameter to rank in Series, DataFrame (:issue:`875`)
- Add sort_columns parameter to allow unsorted plots (:issue:`918`)
- IPython tab completion on GroupBy objects

API Changes
~~~~~~~~~~~

- Series.sum returns 0 instead of NA when called on an empty
  series. Analogously for a DataFrame whose rows or columns are length 0
  (:issue:`844`)

Improvements to existing features
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- Don't use groups dict in Grouper.size (:issue:`860`)
- Use khash for Series.value_counts, add raw function to algorithms.py (:issue:`861`)
- Enable column access via attributes on GroupBy (:issue:`882`)
- Enable setting existing columns (only) via attributes on DataFrame, Panel
  (:issue:`883`)
- Intercept __builtin__.sum in groupby (:issue:`885`)
- Can pass dict to DataFrame.fillna to use different values per column (:issue:`661`)
- Can select multiple hierarchical groups by passing list of values in .ix
  (:issue:`134`)
- Add level keyword to ``drop`` for dropping values from a level (:issue:`159`)
- Add ``coerce_float`` option on DataFrame.from_records (:issue:`893`)
- Raise exception if passed date_parser fails in ``read_csv``
- Add ``axis`` option to DataFrame.fillna (:issue:`174`)
- Fixes to Panel to make it easier to subclass (:issue:`888`)

Bug Fixes
~~~~~~~~~

- Fix overflow-related bugs in groupby (:issue:`850`, :issue:`851`)
- Fix unhelpful error message in parsers (:issue:`856`)
- Better err msg for failed boolean slicing of dataframe (:issue:`859`)
- Series.count cannot accept a string (level name) in the level argument (:issue:`869`)
- Group index platform int check (:issue:`870`)
- concat on axis=1 and ignore_index=True raises TypeError (:issue:`871`)
- Further unicode handling issues resolved (:issue:`795`)
- Fix failure in multiindex-based access in Panel (:issue:`880`)
- Fix DataFrame boolean slice assignment failure (:issue:`881`)
- Fix combineAdd NotImplementedError for SparseDataFrame (:issue:`887`)
- Fix DataFrame.to_html encoding and columns (:issue:`890`, :issue:`891`, :issue:`909`)
- Fix na-filling handling in mixed-type DataFrame (:issue:`910`)
- Fix to DataFrame.set_value with non-existent row/col (:issue:`911`)
- Fix malformed block in groupby when excluding nuisance columns (:issue:`916`)
- Fix inconsistent NA handling in dtype=object arrays (:issue:`925`)
- Fix missing center-of-mass computation in ewmcov (:issue:`862`)
- Don't raise exception when opening read-only HDF5 file (:issue:`847`)
- Fix possible out-of-bounds memory access in 0-length Series (:issue:`917`)

pandas 0.7.1
------------

**Release date:** February 29, 2012

New Features
~~~~~~~~~~~~

- Add ``to_clipboard`` function to pandas namespace for writing objects to
  the system clipboard (:issue:`774`)
- Add ``itertuples`` method to DataFrame for iterating through the rows of a
  dataframe as tuples (:issue:`818`)
- Add ability to pass fill_value and method to DataFrame and Series align
  method (:issue:`806`, :issue:`807`)
- Add fill_value option to reindex, align methods (:issue:`784`)
- Enable concat to produce DataFrame from Series (:issue:`787`)
- Add ``between`` method to Series (:issue:`802`)
- Add HTML representation hook to DataFrame for the IPython HTML notebook
  (:issue:`773`)
- Support for reading Excel 2007 XML documents using openpyxl

Improvements to existing features
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- Improve performance and memory usage of fillna on DataFrame
- Can concatenate a list of Series along axis=1 to obtain a DataFrame (:issue:`787`)

Bug Fixes
~~~~~~~~~

- Fix memory leak when inserting large number of columns into a single
  DataFrame (:issue:`790`)
- Appending length-0 DataFrame with new columns would not result in those new
  columns being part of the resulting concatenated DataFrame (:issue:`782`)
- Fixed groupby corner case when passing dictionary grouper and as_index is
  False (:issue:`819`)
- Fixed bug whereby bool array sometimes had object dtype (:issue:`820`)
- Fix exception thrown on np.diff (:issue:`816`)
- Fix to_records where columns are non-strings (:issue:`822`)
- Fix Index.intersection where indices have incomparable types (:issue:`811`)
- Fix ExcelFile throwing an exception for two-line file (:issue:`837`)
- Add clearer error message in csv parser (:issue:`835`)
- Fix loss of fractional seconds in HDFStore (:issue:`513`)
- Fix DataFrame join where columns have datetimes (:issue:`787`)
- Work around NumPy performance issue in take (:issue:`817`)
- Improve comparison operations for NA-friendliness (:issue:`801`)
- Fix indexing operation for floating point values (:issue:`780`, :issue:`798`)
- Fix groupby case resulting in malformed dataframe (:issue:`814`)
- Fix behavior of reindex of Series dropping name (:issue:`812`)
- Improve on redundant groupby computation (:issue:`775`)
- Catch possible NA assignment to int/bool series with exception (:issue:`839`)

pandas 0.7.0
------------

**Release date:** 2/9/2012

New Features
~~~~~~~~~~~~

- New ``merge`` function for efficiently performing full gamut of database /
  relational-algebra operations. Refactored existing join methods to use the
  new infrastructure, resulting in substantial performance gains (:issue:`220`,
  :issue:`249`, :issue:`267`)
- New ``concat`` function for concatenating DataFrame or Panel objects along
  an axis. Can form union or intersection of the other axes. Improves
  performance of ``DataFrame.append`` (:issue:`468`, :issue:`479`, :issue:`273`)
- Handle differently-indexed output values in ``DataFrame.apply`` (:issue:`498`)
- Can pass list of dicts (e.g., a list of shallow JSON objects) to DataFrame
  constructor (:issue:`526`)
- Add ``reorder_levels`` method to Series and DataFrame (:issue:`534`)
- Add dict-like ``get`` function to DataFrame and Panel (:issue:`521`)
- ``DataFrame.iterrows`` method for efficiently iterating through the rows of
  a DataFrame
- Added ``DataFrame.to_panel`` with code adapted from ``LongPanel.to_long``
- ``reindex_axis`` method added to DataFrame
- Add ``level`` option to binary arithmetic functions on ``DataFrame`` and
  ``Series``
- Add ``level`` option to the ``reindex`` and ``align`` methods on Series and
  DataFrame for broadcasting values across a level (:issue:`542`, :issue:`552`, others)
- Add attribute-based item access to ``Panel`` and add IPython completion (PR
  :issue:`554`)
- Add ``logy`` option to ``Series.plot`` for log-scaling on the Y axis
- Add ``index``, ``header``, and ``justify`` options to
  ``DataFrame.to_string``. Add option to   (:issue:`570`, :issue:`571`)
- Can pass multiple DataFrames to ``DataFrame.join`` to join on index (:issue:`115`)
- Can pass multiple Panels to ``Panel.join`` (:issue:`115`)
- Can pass multiple DataFrames to `DataFrame.append` to concatenate (stack)
  and multiple Series to ``Series.append`` too
- Added ``justify`` argument to ``DataFrame.to_string`` to allow different
  alignment of column headers
- Add ``sort`` option to GroupBy to allow disabling sorting of the group keys
  for potential speedups (:issue:`595`)
- Can pass MaskedArray to Series constructor (:issue:`563`)
- Add Panel item access via attributes and IPython completion (:issue:`554`)
- Implement ``DataFrame.lookup``, fancy-indexing analogue for retrieving
  values given a sequence of row and column labels (:issue:`338`)
- Add ``verbose`` option to ``read_csv`` and ``read_table`` to show number of
  NA values inserted in non-numeric columns (:issue:`614`)
- Can pass a list of dicts or Series to ``DataFrame.append`` to concatenate
  multiple rows (:issue:`464`)
- Add ``level`` argument to ``DataFrame.xs`` for selecting data from other
  MultiIndex levels. Can take one or more levels with potentially a tuple of
  keys for flexible retrieval of data (:issue:`371`, :issue:`629`)
- New ``crosstab`` function for easily computing frequency tables (:issue:`170`)
- Can pass a list of functions to aggregate with groupby on a DataFrame,
  yielding an aggregated result with hierarchical columns (:issue:`166`)
- Add integer-indexing functions ``iget`` in Series and ``irow`` / ``iget``
  in DataFrame (:issue:`628`)
- Add new ``Series.unique`` function, significantly faster than
  ``numpy.unique`` (:issue:`658`)
- Add new ``cummin`` and ``cummax`` instance methods to ``Series`` and
  ``DataFrame`` (:issue:`647`)
- Add new ``value_range`` function to return min/max of a dataframe (:issue:`288`)
- Add ``drop`` parameter to ``reset_index`` method of ``DataFrame`` and added
  method to ``Series`` as well (:issue:`699`)
- Add ``isin`` method to Index objects, works just like ``Series.isin`` (GH
  :issue:`657`)
- Implement array interface on Panel so that ufuncs work (re: :issue:`740`)
- Add ``sort`` option to ``DataFrame.join`` (:issue:`731`)
- Improved handling of NAs (propagation) in binary operations with
  dtype=object arrays (:issue:`737`)
- Add ``abs`` method to Pandas objects
- Added ``algorithms`` module to start collecting central algos

API Changes
~~~~~~~~~~~

- Label-indexing with integer indexes now raises KeyError if a label is not
  found instead of falling back on location-based indexing (:issue:`700`)
- Label-based slicing via ``ix`` or ``[]`` on Series will now only work if
  exact matches for the labels are found or if the index is monotonic (for
  range selections)
- Label-based slicing and sequences of labels can be passed to ``[]`` on a
  Series for both getting and setting (:issue:`86`)
- `[]` operator (``__getitem__`` and ``__setitem__``) will raise KeyError
  with integer indexes when an index is not contained in the index. The prior
  behavior would fall back on position-based indexing if a key was not found
  in the index which would lead to subtle bugs. This is now consistent with
  the behavior of ``.ix`` on DataFrame and friends (:issue:`328`)
- Rename ``DataFrame.delevel`` to ``DataFrame.reset_index`` and add
  deprecation warning
- `Series.sort` (an in-place operation) called on a Series which is a view on
  a larger array (e.g. a column in a DataFrame) will generate an Exception to
  prevent accidentally modifying the data source (:issue:`316`)
- Refactor to remove deprecated ``LongPanel`` class (:issue:`552`)
- Deprecated ``Panel.to_long``, renamed to ``to_frame``
- Deprecated ``colSpace`` argument in ``DataFrame.to_string``, renamed to
  ``col_space``
- Rename ``precision`` to ``accuracy`` in engineering float formatter (GH
  :issue:`395`)
- The default delimiter for ``read_csv`` is comma rather than letting
  ``csv.Sniffer`` infer it
- Rename ``col_or_columns`` argument in ``DataFrame.drop_duplicates`` (GH
  :issue:`734`)

Improvements to existing features
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- Better error message in DataFrame constructor when passed column labels
  don't match data (:issue:`497`)
- Substantially improve performance of multi-GroupBy aggregation when a
  Python function is passed, reuse ndarray object in Cython (:issue:`496`)
- Can store objects indexed by tuples and floats in HDFStore (:issue:`492`)
- Don't print length by default in Series.to_string, add `length` option (GH
  :issue:`489`)
- Improve Cython code for multi-groupby to aggregate without having to sort
  the data (:issue:`93`)
- Improve MultiIndex reindexing speed by storing tuples in the MultiIndex,
  test for backwards unpickling compatibility
- Improve column reindexing performance by using specialized Cython take
  function
- Further performance tweaking of Series.__getitem__ for standard use cases
- Avoid Index dict creation in some cases (i.e. when getting slices, etc.),
  regression from prior versions
- Friendlier error message in setup.py if NumPy not installed
- Use common set of NA-handling operations (sum, mean, etc.) in Panel class
  also (:issue:`536`)
- Default name assignment when calling ``reset_index`` on DataFrame with a
  regular (non-hierarchical) index (:issue:`476`)
- Use Cythonized groupers when possible in Series/DataFrame stat ops with
  ``level`` parameter passed (:issue:`545`)
- Ported skiplist data structure to C to speed up ``rolling_median`` by about
  5-10x in most typical use cases (:issue:`374`)
- Some performance enhancements in constructing a Panel from a dict of
  DataFrame objects
- Made ``Index._get_duplicates`` a public method by removing the underscore
- Prettier printing of floats, and column spacing fix (:issue:`395`, :issue:`571`)
- Add ``bold_rows`` option to DataFrame.to_html (:issue:`586`)
- Improve the performance of ``DataFrame.sort_index`` by up to 5x or more
  when sorting by multiple columns
- Substantially improve performance of DataFrame and Series constructors when
  passed a nested dict or dict, respectively (:issue:`540`, :issue:`621`)
- Modified setup.py so that pip / setuptools will install dependencies (GH
  :issue:`507`, various pull requests)
- Unstack called on DataFrame with non-MultiIndex will return Series (GH
  :issue:`477`)
- Improve DataFrame.to_string and console formatting to be more consistent in
  the number of displayed digits (:issue:`395`)
- Use bottleneck if available for performing NaN-friendly statistical
  operations that it implemented (:issue:`91`)
- Monkey-patch context to traceback in ``DataFrame.apply`` to indicate which
  row/column the function application failed on (:issue:`614`)
- Improved ability of read_table and read_clipboard to parse
  console-formatted DataFrames (can read the row of index names, etc.)
- Can pass list of group labels (without having to convert to an ndarray
  yourself) to ``groupby`` in some cases (:issue:`659`)
- Use ``kind`` argument to Series.order for selecting different sort kinds
  (:issue:`668`)
- Add option to Series.to_csv to omit the index (:issue:`684`)
- Add ``delimiter`` as an alternative to ``sep`` in ``read_csv`` and other
  parsing functions
- Substantially improved performance of groupby on DataFrames with many
  columns by aggregating blocks of columns all at once (:issue:`745`)
- Can pass a file handle or StringIO to Series/DataFrame.to_csv (:issue:`765`)
- Can pass sequence of integers to DataFrame.irow(icol) and Series.iget, (GH
  :issue:`654`)
- Prototypes for some vectorized string functions
- Add float64 hash table to solve the Series.unique problem with NAs (:issue:`714`)
- Memoize objects when reading from file to reduce memory footprint
- Can get and set a column of a DataFrame with hierarchical columns
  containing "empty" ('') lower levels without passing the empty levels (PR
  :issue:`768`)

Bug Fixes
~~~~~~~~~

- Raise exception in out-of-bounds indexing of Series instead of
  seg-faulting, regression from earlier releases (:issue:`495`)
- Fix error when joining DataFrames of different dtypes within the same
  type class (e.g. float32 and float64) (:issue:`486`)
- Fix bug in Series.min/Series.max on objects like datetime.datetime (GH
  :issue:`487`)
- Preserve index names in Index.union (:issue:`501`)
- Fix bug in Index joining causing subclass information (like DateRange type)
  to be lost in some cases (:issue:`500`)
- Accept empty list as input to DataFrame constructor, regression from 0.6.0
  (:issue:`491`)
- Can output DataFrame and Series with ndarray objects in a dtype=object
  array (:issue:`490`)
- Return empty string from Series.to_string when called on empty Series (GH
  :issue:`488`)
- Fix exception passing empty list to DataFrame.from_records
- Fix Index.format bug (excluding name field) with datetimes with time info
- Fix scalar value access in Series to always return NumPy scalars,
  regression from prior versions (:issue:`510`)
- Handle rows skipped at beginning of file in read_* functions (:issue:`505`)
- Handle improper dtype casting in ``set_value`` methods
- Unary '-' / __neg__ operator on DataFrame was returning integer values
- Unbox 0-dim ndarrays from certain operators like all, any in Series
- Fix handling of missing columns (was combine_first-specific) in
  DataFrame.combine for general case (:issue:`529`)
- Fix type inference logic with boolean lists and arrays in DataFrame indexing
- Use centered sum of squares in R-square computation if entity_effects=True
  in panel regression
- Handle all NA case in Series.{corr, cov}, was raising exception (:issue:`548`)
- Aggregating by multiple levels with ``level`` argument to DataFrame, Series
  stat method, was broken (:issue:`545`)
- Fix Cython buf when converter passed to read_csv produced a numeric array
  (buffer dtype mismatch when passed to Cython type inference function) (GH
  :issue:`546`)
- Fix exception when setting scalar value using .ix on a DataFrame with a
  MultiIndex (:issue:`551`)
- Fix outer join between two DateRanges with different offsets that returned
  an invalid DateRange
- Cleanup DataFrame.from_records failure where index argument is an integer
- Fix Data.from_records failure when passed a dictionary
- Fix NA handling in {Series, DataFrame}.rank with non-floating point dtypes
- Fix bug related to integer type-checking in .ix-based indexing
- Handle non-string index name passed to DataFrame.from_records
- DataFrame.insert caused the columns name(s) field to be discarded (:issue:`527`)
- Fix erroneous in monotonic many-to-one left joins
- Fix DataFrame.to_string to remove extra column white space (:issue:`571`)
- Format floats to default to same number of digits (:issue:`395`)
- Added decorator to copy docstring from one function to another (:issue:`449`)
- Fix error in monotonic many-to-one left joins
- Fix __eq__ comparison between DateOffsets with different relative delta
  keywords passed
- Fix exception caused by parser converter returning strings (:issue:`583`)
- Fix MultiIndex formatting bug with integer names (:issue:`601`)
- Fix bug in handling of non-numeric aggregates in Series.groupby (:issue:`612`)
- Fix TypeError with tuple subclasses (e.g. namedtuple) in
  DataFrame.from_records (:issue:`611`)
- Catch misreported console size when running IPython within Emacs
- Fix minor bug in pivot table margins, loss of index names and length-1
  'All' tuple in row labels
- Add support for legacy WidePanel objects to be read from HDFStore
- Fix out-of-bounds segfault in pad_object and backfill_object methods when
  either source or target array are empty
- Could not create a new column in a DataFrame from a list of tuples
- Fix bugs preventing SparseDataFrame and SparseSeries working with groupby
  (:issue:`666`)
- Use sort kind in Series.sort / argsort (:issue:`668`)
- Fix DataFrame operations on non-scalar, non-pandas objects (:issue:`672`)
- Don't convert DataFrame column to integer type when passing integer to
  __setitem__ (:issue:`669`)
- Fix downstream bug in pivot_table caused by integer level names in
  MultiIndex (:issue:`678`)
- Fix SparseSeries.combine_first when passed a dense Series (:issue:`687`)
- Fix performance regression in HDFStore loading when DataFrame or Panel
  stored in table format with datetimes
- Raise Exception in DateRange when offset with n=0 is passed (:issue:`683`)
- Fix get/set inconsistency with .ix property and integer location but
  non-integer index (:issue:`707`)
- Use right dropna function for SparseSeries. Return dense Series for NA fill
  value (:issue:`730`)
- Fix Index.format bug causing incorrectly string-formatted Series with
  datetime indexes (:issue:`726`, :issue:`758`)
- Fix errors caused by object dtype arrays passed to ols (:issue:`759`)
- Fix error where column names lost when passing list of labels to
  DataFrame.__getitem__, (:issue:`662`)
- Fix error whereby top-level week iterator overwrote week instance
- Fix circular reference causing memory leak in sparse array / series /
  frame, (:issue:`663`)
- Fix integer-slicing from integers-as-floats (:issue:`670`)
- Fix zero division errors in nanops from object dtype arrays in all NA case
  (:issue:`676`)
- Fix csv encoding when using unicode (:issue:`705`, :issue:`717`, :issue:`738`)
- Fix assumption that each object contains every unique block type in concat,
  (:issue:`708`)
- Fix sortedness check of multiindex in to_panel (:issue:`719`, 720)
- Fix that None was not treated as NA in PyObjectHashtable
- Fix hashing dtype because of endianness confusion (:issue:`747`, :issue:`748`)
- Fix SparseSeries.dropna to return dense Series in case of NA fill value (GH
  :issue:`730`)
- Use map_infer instead of np.vectorize. handle NA sentinels if converter
  yields numeric array, (:issue:`753`)
- Fixes and improvements to DataFrame.rank (:issue:`742`)
- Fix catching AttributeError instead of NameError for bottleneck
- Try to cast non-MultiIndex to better dtype when calling reset_index (:issue:`726`
  :issue:`440`)
- Fix #1.QNAN0' float bug on 2.6/win64
- Allow subclasses of dicts in DataFrame constructor, with tests
- Fix problem whereby set_index destroys column multiindex (:issue:`764`)
- Hack around bug in generating DateRange from naive DateOffset (:issue:`770`)
- Fix bug in DateRange.intersection causing incorrect results with some
  overlapping ranges (:issue:`771`)

Thanks
~~~~~~

- Craig Austin
- Chris Billington
- Marius Cobzarenco
- Mario Gamboa-Cavazos
- Hans-Martin Gaudecker
- Arthur Gerigk
- Yaroslav Halchenko
- Jeff Hammerbacher
- Matt Harrison
- Andreas Hilboll
- Luc Kesters
- Adam Klein
- Gregg Lind
- Solomon Negusse
- Wouter Overmeire
- Christian Prinoth
- Jeff Reback
- Sam Reckoner
- Craig Reeson
- Jan Schulz
- Skipper Seabold
- Ted Square
- Graham Taylor
- Aman Thakral
- Chris Uga
- Dieter Vandenbussche
- Texas P.
- Pinxing Ye
- ... and everyone I forgot

pandas 0.6.1
------------

**Release date:** 12/13/2011

API Changes
~~~~~~~~~~~

- Rename `names` argument in DataFrame.from_records to `columns`. Add
  deprecation warning
- Boolean get/set operations on Series with boolean Series will reindex
  instead of requiring that the indexes be exactly equal (:issue:`429`)

New Features
~~~~~~~~~~~~

- Can pass Series to DataFrame.append with ignore_index=True for appending a
  single row (:issue:`430`)
- Add Spearman and Kendall correlation options to Series.corr and
  DataFrame.corr (:issue:`428`)
- Add new `get_value` and `set_value` methods to Series, DataFrame, and Panel
  to very low-overhead access to scalar elements. df.get_value(row, column)
  is about 3x faster than df[column][row] by handling fewer cases (:issue:`437`,
  :issue:`438`). Add similar methods to sparse data structures for compatibility
- Add Qt table widget to sandbox (:issue:`435`)
- DataFrame.align can accept Series arguments, add axis keyword (:issue:`461`)
- Implement new SparseList and SparseArray data structures. SparseSeries now
  derives from SparseArray (:issue:`463`)
- max_columns / max_rows options in set_printoptions (:issue:`453`)
- Implement Series.rank and DataFrame.rank, fast versions of
  scipy.stats.rankdata (:issue:`428`)
- Implement DataFrame.from_items alternate constructor (:issue:`444`)
- DataFrame.convert_objects method for inferring better dtypes for object
  columns (:issue:`302`)
- Add rolling_corr_pairwise function for computing Panel of correlation
  matrices (:issue:`189`)
- Add `margins` option to `pivot_table` for computing subgroup aggregates (GH
  :issue:`114`)
- Add `Series.from_csv` function (:issue:`482`)

Improvements to existing features
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- Improve memory usage of `DataFrame.describe` (do not copy data
  unnecessarily) (:issue:`425`)
- Use same formatting function for outputting floating point Series to console
  as in DataFrame (:issue:`420`)
- DataFrame.delevel will try to infer better dtype for new columns (:issue:`440`)
- Exclude non-numeric types in DataFrame.{corr, cov}
- Override Index.astype to enable dtype casting (:issue:`412`)
- Use same float formatting function for Series.__repr__ (:issue:`420`)
- Use available console width to output DataFrame columns (:issue:`453`)
- Accept ndarrays when setting items in Panel (:issue:`452`)
- Infer console width when printing __repr__ of DataFrame to console (PR
  :issue:`453`)
- Optimize scalar value lookups in the general case by 25% or more in Series
  and DataFrame
- Can pass DataFrame/DataFrame and DataFrame/Series to
  rolling_corr/rolling_cov (:issue:`462`)
- Fix performance regression in cross-sectional count in DataFrame, affecting
  DataFrame.dropna speed
- Column deletion in DataFrame copies no data (computes views on blocks) (GH
  :issue:`158`)
- MultiIndex.get_level_values can take the level name
- More helpful error message when DataFrame.plot fails on one of the columns
  (:issue:`478`)
- Improve performance of DataFrame.{index, columns} attribute lookup

Bug Fixes
~~~~~~~~~

- Fix O(K^2) memory leak caused by inserting many columns without
  consolidating, had been present since 0.4.0 (:issue:`467`)
- `DataFrame.count` should return Series with zero instead of NA with length-0
  axis (:issue:`423`)
- Fix Yahoo! Finance API usage in pandas.io.data (:issue:`419`, :issue:`427`)
- Fix upstream bug causing failure in Series.align with empty Series (:issue:`434`)
- Function passed to DataFrame.apply can return a list, as long as it's the
  right length. Regression from 0.4 (:issue:`432`)
- Don't "accidentally" upcast scalar values when indexing using .ix (:issue:`431`)
- Fix groupby exception raised with as_index=False and single column selected
  (:issue:`421`)
- Implement DateOffset.__ne__ causing downstream bug (:issue:`456`)
- Fix __doc__-related issue when converting py -> pyo with py2exe
- Bug fix in left join Cython code with duplicate monotonic labels
- Fix bug when unstacking multiple levels described in :issue:`451`
- Exclude NA values in dtype=object arrays, regression from 0.5.0 (:issue:`469`)
- Use Cython map_infer function in DataFrame.applymap to properly infer
  output type, handle tuple return values and other things that were breaking
  (:issue:`465`)
- Handle floating point index values in HDFStore (:issue:`454`)
- Fixed stale column reference bug (cached Series object) caused by type
  change / item deletion in DataFrame (:issue:`473`)
- Index.get_loc should always raise Exception when there are duplicates
- Handle differently-indexed Series input to DataFrame constructor (:issue:`475`)
- Omit nuisance columns in multi-groupby with Python function
- Buglet in handling of single grouping in general apply
- Handle type inference properly when passing list of lists or tuples to
  DataFrame constructor (:issue:`484`)
- Preserve Index / MultiIndex names in GroupBy.apply concatenation step (GH
  :issue:`481`)

Thanks
~~~~~~

- Ralph Bean
- Luca Beltrame
- Marius Cobzarenco
- Andreas Hilboll
- Jev Kuznetsov
- Adam Lichtenstein
- Wouter Overmeire
- Fernando Perez
- Nathan Pinger
- Christian Prinoth
- Alex Reyfman
- Joon Ro
- Chang She
- Ted Square
- Chris Uga
- Dieter Vandenbussche

pandas 0.6.0
------------

**Release date:** 11/25/2011

API Changes
~~~~~~~~~~~

- Arithmetic methods like `sum` will attempt to sum dtype=object values by
  default instead of excluding them (:issue:`382`)

New Features
~~~~~~~~~~~~

- Add `melt` function to `pandas.core.reshape`
- Add `level` parameter to group by level in Series and DataFrame
  descriptive statistics (:issue:`313`)
- Add `head` and `tail` methods to Series, analogous to DataFrame (PR
  :issue:`296`)
- Add `Series.isin` function which checks if each value is contained in a
  passed sequence (:issue:`289`)
- Add `float_format` option to `Series.to_string`
- Add `skip_footer` (:issue:`291`) and `converters` (:issue:`343`) options to
  `read_csv` and `read_table`
- Add proper, tested weighted least squares to standard and panel OLS (GH
  :issue:`303`)
- Add `drop_duplicates` and `duplicated` functions for removing duplicate
  DataFrame rows and checking for duplicate rows, respectively (:issue:`319`)
- Implement logical (boolean) operators ``&``, ``|``, ``^`` on DataFrame
  (:issue:`347`)
- Add `Series.mad`, mean absolute deviation, matching DataFrame
- Add `QuarterEnd` DateOffset (:issue:`321`)
- Add matrix multiplication function `dot` to DataFrame (:issue:`65`)
- Add `orient` option to `Panel.from_dict` to ease creation of mixed-type
  Panels (:issue:`359`, :issue:`301`)
- Add `DataFrame.from_dict` with similar `orient` option
- Can now pass list of tuples or list of lists to `DataFrame.from_records`
  for fast conversion to DataFrame (:issue:`357`)
- Can pass multiple levels to groupby, e.g. `df.groupby(level=[0, 1])` (GH
  :issue:`103`)
- Can sort by multiple columns in `DataFrame.sort_index` (:issue:`92`, :issue:`362`)
- Add fast `get_value` and `put_value` methods to DataFrame and
  micro-performance tweaks (:issue:`360`)
- Add `cov` instance methods to Series and DataFrame (:issue:`194`, :issue:`362`)
- Add bar plot option to `DataFrame.plot` (:issue:`348`)
- Add `idxmin` and `idxmax` functions to Series and DataFrame for computing
  index labels achieving maximum and minimum values (:issue:`286`)
- Add `read_clipboard` function for parsing DataFrame from OS clipboard,
  should work across platforms (:issue:`300`)
- Add `nunique` function to Series for counting unique elements (:issue:`297`)
- DataFrame constructor will use Series name if no columns passed (:issue:`373`)
- Support regular expressions and longer delimiters in read_table/read_csv,
  but does not handle quoted strings yet (:issue:`364`)
- Add `DataFrame.to_html` for formatting DataFrame to HTML (:issue:`387`)
- MaskedArray can be passed to DataFrame constructor and masked values will be
  converted to NaN (:issue:`396`)
- Add `DataFrame.boxplot` function (:issue:`368`, others)
- Can pass extra args, kwds to DataFrame.apply (:issue:`376`)

Improvements to existing features
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- Raise more helpful exception if date parsing fails in DateRange (:issue:`298`)
- Vastly improved performance of GroupBy on axes with a MultiIndex (:issue:`299`)
- Print level names in hierarchical index in Series repr (:issue:`305`)
- Return DataFrame when performing GroupBy on selected column and
  as_index=False (:issue:`308`)
- Can pass vector to `on` argument in `DataFrame.join` (:issue:`312`)
- Don't show Series name if it's None in the repr, also omit length for short
  Series (:issue:`317`)
- Show legend by default in `DataFrame.plot`, add `legend` boolean flag (GH
  :issue:`324`)
- Significantly improved performance of `Series.order`, which also makes
  np.unique called on a Series faster (:issue:`327`)
- Faster cythonized count by level in Series and DataFrame (:issue:`341`)
- Raise exception if dateutil 2.0 installed on Python 2.x runtime (:issue:`346`)
- Significant GroupBy performance enhancement with multiple keys with many
  "empty" combinations
- New Cython vectorized function `map_infer` speeds up `Series.apply` and
  `Series.map` significantly when passed elementwise Python function,
  motivated by :issue:`355`
- Cythonized `cache_readonly`, resulting in substantial micro-performance
  enhancements throughout the code base (:issue:`361`)
- Special Cython matrix iterator for applying arbitrary reduction operations
  with 3-5x better performance than `np.apply_along_axis` (:issue:`309`)
- Add `raw` option to `DataFrame.apply` for getting better performance when
  the passed function only requires an ndarray (:issue:`309`)
- Improve performance of `MultiIndex.from_tuples`
- Can pass multiple levels to `stack` and `unstack` (:issue:`370`)
- Can pass multiple values columns to `pivot_table` (:issue:`381`)
- Can call `DataFrame.delevel` with standard Index with name set (:issue:`393`)
- Use Series name in GroupBy for result index (:issue:`363`)
- Refactor Series/DataFrame stat methods to use common set of NaN-friendly
  function
- Handle NumPy scalar integers at C level in Cython conversion routines

Bug Fixes
~~~~~~~~~

- Fix bug in `DataFrame.to_csv` when writing a DataFrame with an index
  name (:issue:`290`)
- DataFrame should clear its Series caches on consolidation, was causing
  "stale" Series to be returned in some corner cases (:issue:`304`)
- DataFrame constructor failed if a column had a list of tuples (:issue:`293`)
- Ensure that `Series.apply` always returns a Series and implement
  `Series.round` (:issue:`314`)
- Support boolean columns in Cythonized groupby functions (:issue:`315`)
- `DataFrame.describe` should not fail if there are no numeric columns,
  instead return categorical describe (:issue:`323`)
- Fixed bug which could cause columns to be printed in wrong order in
  `DataFrame.to_string` if specific list of columns passed (:issue:`325`)
- Fix legend plotting failure if DataFrame columns are integers (:issue:`326`)
- Shift start date back by one month for Yahoo! Finance API in pandas.io.data
  (:issue:`329`)
- Fix `DataFrame.join` failure on unconsolidated inputs (:issue:`331`)
- DataFrame.min/max will no longer fail on mixed-type DataFrame (:issue:`337`)
- Fix `read_csv` / `read_table` failure when passing list to index_col that is
  not in ascending order (:issue:`349`)
- Fix failure passing Int64Index to Index.union when both are monotonic
- Fix error when passing SparseSeries to (dense) DataFrame constructor
- Added missing bang at top of setup.py (:issue:`352`)
- Change `is_monotonic` on MultiIndex so it properly compares the tuples
- Fix MultiIndex outer join logic (:issue:`351`)
- Set index name attribute with single-key groupby (:issue:`358`)
- Bug fix in reflexive binary addition in Series and DataFrame for
  non-commutative operations (like string concatenation) (:issue:`353`)
- setupegg.py will invoke Cython (:issue:`192`)
- Fix block consolidation bug after inserting column into MultiIndex (:issue:`366`)
- Fix bug in join operations between Index and Int64Index (:issue:`367`)
- Handle min_periods=0 case in moving window functions (:issue:`365`)
- Fixed corner cases in DataFrame.apply/pivot with empty DataFrame (:issue:`378`)
- Fixed repr exception when Series name is a tuple
- Always return DateRange from `asfreq` (:issue:`390`)
- Pass level names to `swaplavel` (:issue:`379`)
- Don't lose index names in `MultiIndex.droplevel` (:issue:`394`)
- Infer more proper return type in `DataFrame.apply` when no columns or rows
  depending on whether the passed function is a reduction (:issue:`389`)
- Always return NA/NaN from Series.min/max and DataFrame.min/max when all of a
  row/column/values are NA (:issue:`384`)
- Enable partial setting with .ix / advanced indexing (:issue:`397`)
- Handle mixed-type DataFrames correctly in unstack, do not lose type
  information (:issue:`403`)
- Fix integer name formatting bug in Index.format and in Series.__repr__
- Handle label types other than string passed to groupby (:issue:`405`)
- Fix bug in .ix-based indexing with partial retrieval when a label is not
  contained in a level
- Index name was not being pickled (:issue:`408`)
- Level name should be passed to result index in GroupBy.apply (:issue:`416`)

Thanks
~~~~~~

- Craig Austin
- Marius Cobzarenco
- Joel Cross
- Jeff Hammerbacher
- Adam Klein
- Thomas Kluyver
- Jev Kuznetsov
- Kieran O'Mahony
- Wouter Overmeire
- Nathan Pinger
- Christian Prinoth
- Skipper Seabold
- Chang She
- Ted Square
- Aman Thakral
- Chris Uga
- Dieter Vandenbussche
- carljv
- rsamson

pandas 0.5.0
------------

**Release date:** 10/24/2011

This release of pandas includes a number of API changes (see below) and cleanup of deprecated APIs
from pre-0.4.0 releases. There are also bug fixes, new features, numerous significant performance enhancements, and includes a new ipython
completer hook to enable tab completion of DataFrame columns accesses and attributes (a new feature).

In addition to the changes listed here from 0.4.3 to 0.5.0, the minor releases 4.1,
0.4.2, and 0.4.3 brought some significant new functionality and performance improvements that are worth taking a look at.

Thanks to all for bug reports, contributed patches and generally providing feedback on the library.

API Changes
~~~~~~~~~~~

- `read_table`, `read_csv`, and `ExcelFile.parse` default arguments for
  `index_col` is now None. To use one or more of the columns as the resulting
  DataFrame's index, these must be explicitly specified now
- Parsing functions like `read_csv` no longer parse dates by default (GH
  :issue:`225`)
- Removed `weights` option in panel regression which was not doing anything
  principled (:issue:`155`)
- Changed `buffer` argument name in `Series.to_string` to `buf`
- `Series.to_string` and `DataFrame.to_string` now return strings by default
  instead of printing to sys.stdout
- Deprecated `nanRep` argument in various `to_string` and `to_csv` functions
  in favor of `na_rep`. Will be removed in 0.6 (:issue:`275`)
- Renamed `delimiter` to `sep` in `DataFrame.from_csv` for consistency
- Changed order of `Series.clip` arguments to match those of `numpy.clip` and
  added (unimplemented) `out` argument so `numpy.clip` can be called on a
  Series (:issue:`272`)
- Series functions renamed (and thus deprecated) in 0.4 series have been
  removed:

  - `asOf`, use `asof`
  - `toDict`, use `to_dict`
  - `toString`, use `to_string`
  - `toCSV`, use `to_csv`
  - `merge`, use `map`
  - `applymap`, use `apply`
  - `combineFirst`, use `combine_first`
  - `_firstTimeWithValue` use `first_valid_index`
  - `_lastTimeWithValue` use `last_valid_index`

- DataFrame functions renamed / deprecated in 0.4 series have been removed:

  - `asMatrix` method, use `as_matrix` or `values` attribute
  - `combineFirst`, use `combine_first`
  - `getXS`, use `xs`
  - `merge`, use `join`
  - `fromRecords`, use `from_records`
  - `fromcsv`, use `from_csv`
  - `toRecords`, use `to_records`
  - `toDict`, use `to_dict`
  - `toString`, use `to_string`
  - `toCSV`, use `to_csv`
  - `_firstTimeWithValue` use `first_valid_index`
  - `_lastTimeWithValue` use `last_valid_index`
  - `toDataMatrix` is no longer needed
  - `rows()` method, use `index` attribute
  - `cols()` method, use `columns` attribute
  - `dropEmptyRows()`, use `dropna(how='all')`
  - `dropIncompleteRows()`, use `dropna()`
  - `tapply(f)`, use `apply(f, axis=1)`
  - `tgroupby(keyfunc, aggfunc)`, use `groupby` with `axis=1`

Deprecations Removed
~~~~~~~~~~~~~~~~~~~~

  - `indexField` argument in `DataFrame.from_records`
  - `missingAtEnd` argument in `Series.order`. Use `na_last` instead
  - `Series.fromValue` classmethod, use regular `Series` constructor instead
  - Functions `parseCSV`, `parseText`, and `parseExcel` methods in
    `pandas.io.parsers` have been removed
  - `Index.asOfDate` function
  - `Panel.getMinorXS` (use `minor_xs`) and `Panel.getMajorXS` (use
    `major_xs`)
  - `Panel.toWide`, use `Panel.to_wide` instead

New Features
~~~~~~~~~~~~

- Added `DataFrame.align` method with standard join options
- Added `parse_dates` option to `read_csv` and `read_table` methods to
  optionally try to parse dates in the index columns
- Add `nrows`, `chunksize`, and `iterator` arguments to `read_csv` and
  `read_table`. The last two return a new `TextParser` class capable of
  lazily iterating through chunks of a flat file (:issue:`242`)
- Added ability to join on multiple columns in `DataFrame.join` (:issue:`214`)
- Added private `_get_duplicates` function to `Index` for identifying
  duplicate values more easily
- Added column attribute access to DataFrame, e.g. df.A equivalent to df['A']
  if 'A' is a column in the DataFrame (:issue:`213`)
- Added IPython tab completion hook for DataFrame columns. (:issue:`233`, :issue:`230`)
- Implement `Series.describe` for Series containing objects (:issue:`241`)
- Add inner join option to `DataFrame.join` when joining on key(s) (:issue:`248`)
- Can select set of DataFrame columns by passing a list to `__getitem__` (GH
  :issue:`253`)
- Can use & and | to intersection / union Index objects, respectively (GH
  :issue:`261`)
- Added `pivot_table` convenience function to pandas namespace (:issue:`234`)
- Implemented `Panel.rename_axis` function (:issue:`243`)
- DataFrame will show index level names in console output
- Implemented `Panel.take`
- Add `set_eng_float_format` function for setting alternate DataFrame
  floating point string formatting
- Add convenience `set_index` function for creating a DataFrame index from
  its existing columns

Improvements to existing features
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- Major performance improvements in file parsing functions `read_csv` and
  `read_table`
- Added Cython function for converting tuples to ndarray very fast. Speeds up
  many MultiIndex-related operations
- File parsing functions like `read_csv` and `read_table` will explicitly
  check if a parsed index has duplicates and raise a more helpful exception
  rather than deferring the check until later
- Refactored merging / joining code into a tidy class and disabled unnecessary
  computations in the float/object case, thus getting about 10% better
  performance (:issue:`211`)
- Improved speed of `DataFrame.xs` on mixed-type DataFrame objects by about
  5x, regression from 0.3.0 (:issue:`215`)
- With new `DataFrame.align` method, speeding up binary operations between
  differently-indexed DataFrame objects by 10-25%.
- Significantly sped up conversion of nested dict into DataFrame (:issue:`212`)
- Can pass hierarchical index level name to `groupby` instead of the level
  number if desired (:issue:`223`)
- Add support for different delimiters in `DataFrame.to_csv` (:issue:`244`)
- Add more helpful error message when importing pandas post-installation from
  the source directory (:issue:`250`)
- Significantly speed up DataFrame `__repr__` and `count` on large mixed-type
  DataFrame objects
- Better handling of pyx file dependencies in Cython module build (:issue:`271`)

Bug Fixes
~~~~~~~~~

- `read_csv` / `read_table` fixes

  - Be less aggressive about converting float->int in cases of floating point
    representations of integers like 1.0, 2.0, etc.
  - "True"/"False" will not get correctly converted to boolean
  - Index name attribute will get set when specifying an index column
  - Passing column names should force `header=None` (:issue:`257`)
  - Don't modify passed column names when `index_col` is not None
    (:issue:`258`)
  - Can sniff CSV separator in zip file (since seek is not supported, was
    failing before)

- Worked around matplotlib "bug" in which series[:, np.newaxis] fails. Should
  be reported upstream to matplotlib (:issue:`224`)
- DataFrame.iteritems was not returning Series with the name attribute
  set. Also neither was DataFrame._series
- Can store datetime.date objects in HDFStore (:issue:`231`)
- Index and Series names are now stored in HDFStore
- Fixed problem in which data would get upcasted to object dtype in
  GroupBy.apply operations (:issue:`237`)
- Fixed outer join bug with empty DataFrame (:issue:`238`)
- Can create empty Panel (:issue:`239`)
- Fix join on single key when passing list with 1 entry (:issue:`246`)
- Don't raise Exception on plotting DataFrame with an all-NA column (:issue:`251`,
  :issue:`254`)
- Bug min/max errors when called on integer DataFrames (:issue:`241`)
- `DataFrame.iteritems` and `DataFrame._series` not assigning name attribute
- Panel.__repr__ raised exception on length-0 major/minor axes
- `DataFrame.join` on key with empty DataFrame produced incorrect columns
- Implemented `MultiIndex.diff` (:issue:`260`)
- `Int64Index.take` and `MultiIndex.take` lost name field, fix downstream
  issue :issue:`262`
- Can pass list of tuples to `Series` (:issue:`270`)
- Can pass level name to `DataFrame.stack`
- Support set operations between MultiIndex and Index
- Fix many corner cases in MultiIndex set operations
  - Fix MultiIndex-handling bug with GroupBy.apply when returned groups are not
  indexed the same
- Fix corner case bugs in DataFrame.apply
- Setting DataFrame index did not cause Series cache to get cleared
- Various int32 -> int64 platform-specific issues
- Don't be too aggressive converting to integer when parsing file with
  MultiIndex (:issue:`285`)
- Fix bug when slicing Series with negative indices before beginning

Thanks
~~~~~~

- Thomas Kluyver
- Daniel Fortunov
- Aman Thakral
- Luca Beltrame
- Wouter Overmeire

pandas 0.4.3
------------

**Release date:** 10/9/2011

This is largely a bug fix release from 0.4.2 but also includes a handful of new
and enhanced features. Also, pandas can now be installed and used on Python 3
(thanks Thomas Kluyver!).

New Features
~~~~~~~~~~~~

- Python 3 support using 2to3 (:issue:`200`, Thomas Kluyver)
- Add `name` attribute to `Series` and added relevant logic and tests. Name
  now prints as part of `Series.__repr__`
- Add `name` attribute to standard Index so that stacking / unstacking does
  not discard names and so that indexed DataFrame objects can be reliably
  round-tripped to flat files, pickle, HDF5, etc.
- Add `isnull` and `notnull` as instance methods on Series (:issue:`209`, :issue:`203`)

Improvements to existing features
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- Skip xlrd-related unit tests if not installed
- `Index.append` and `MultiIndex.append` can accept a list of Index objects to
  concatenate together
- Altered binary operations on differently-indexed SparseSeries objects to use
  the integer-based (dense) alignment logic which is faster with a larger
  number of blocks (:issue:`205`)
- Refactored `Series.__repr__` to be a bit more clean and consistent

API Changes
~~~~~~~~~~~

- `Series.describe` and `DataFrame.describe` now bring the 25% and 75%
  quartiles instead of the 10% and 90% deciles. The other outputs have not
  changed
- `Series.toString` will print deprecation warning, has been de-camelCased to
  `to_string`

Bug Fixes
~~~~~~~~~

- Fix broken interaction between `Index` and `Int64Index` when calling
  intersection. Implement `Int64Index.intersection`
- `MultiIndex.sortlevel` discarded the level names (:issue:`202`)
- Fix bugs in groupby, join, and append due to improper concatenation of
  `MultiIndex` objects (:issue:`201`)
- Fix regression from 0.4.1, `isnull` and `notnull` ceased to work on other
  kinds of Python scalar objects like `datetime.datetime`
- Raise more helpful exception when attempting to write empty DataFrame or
  LongPanel to `HDFStore` (:issue:`204`)
- Use stdlib csv module to properly escape strings with commas in
  `DataFrame.to_csv` (:issue:`206`, Thomas Kluyver)
- Fix Python ndarray access in Cython code for sparse blocked index integrity
  check
- Fix bug writing Series to CSV in Python 3 (:issue:`209`)
- Miscellaneous Python 3 bug fixes

Thanks
~~~~~~

- Thomas Kluyver
- rsamson

pandas 0.4.2
------------

**Release date:** 10/3/2011

This is a performance optimization release with several bug fixes. The new
Int64Index and new merging / joining Cython code and related Python
infrastructure are the main new additions

New Features
~~~~~~~~~~~~

- Added fast `Int64Index` type with specialized join, union,
  intersection. Will result in significant performance enhancements for
  int64-based time series (e.g. using NumPy's datetime64 one day) and also
  faster operations on DataFrame objects storing record array-like data.
- Refactored `Index` classes to have a `join` method and associated data
  alignment routines throughout the code base to be able to leverage optimized
  joining / merging routines.
- Added `Series.align` method for aligning two series with choice of join
  method
- Wrote faster Cython data alignment / merging routines resulting in
  substantial speed increases
- Added `is_monotonic` property to `Index` classes with associated Cython
  code to evaluate the monotonicity of the `Index` values
- Add method `get_level_values` to `MultiIndex`
- Implemented shallow copy of `BlockManager` object in `DataFrame` internals

Improvements to existing features
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- Improved performance of `isnull` and `notnull`, a regression from v0.3.0
  (:issue:`187`)
- Wrote templating / code generation script to auto-generate Cython code for
  various functions which need to be available for the 4 major data types
  used in pandas (float64, bool, object, int64)
- Refactored code related to `DataFrame.join` so that intermediate aligned
  copies of the data in each `DataFrame` argument do not need to be
  created. Substantial performance increases result (:issue:`176`)
- Substantially improved performance of generic `Index.intersection` and
  `Index.union`
- Improved performance of `DateRange.union` with overlapping ranges and
  non-cacheable offsets (like Minute). Implemented analogous fast
  `DateRange.intersection` for overlapping ranges.
- Implemented `BlockManager.take` resulting in significantly faster `take`
  performance on mixed-type `DataFrame` objects (:issue:`104`)
- Improved performance of `Series.sort_index`
- Significant groupby performance enhancement: removed unnecessary integrity
  checks in DataFrame internals that were slowing down slicing operations to
  retrieve groups
- Added informative Exception when passing dict to DataFrame groupby
  aggregation with axis != 0

API Changes
~~~~~~~~~~~

Bug Fixes
~~~~~~~~~

- Fixed minor unhandled exception in Cython code implementing fast groupby
  aggregation operations
- Fixed bug in unstacking code manifesting with more than 3 hierarchical
  levels
- Throw exception when step specified in label-based slice (:issue:`185`)
- Fix isnull to correctly work with np.float32. Fix upstream bug described in
  :issue:`182`
- Finish implementation of as_index=False in groupby for DataFrame
  aggregation (:issue:`181`)
- Raise SkipTest for pre-epoch HDFStore failure. Real fix will be sorted out
  via datetime64 dtype

Thanks
~~~~~~

- Uri Laserson
- Scott Sinclair

pandas 0.4.1
------------

**Release date:** 9/25/2011

This is primarily a bug fix release but includes some new features and
improvements

New Features
~~~~~~~~~~~~

- Added new `DataFrame` methods `get_dtype_counts` and property `dtypes`
- Setting of values using ``.ix`` indexing attribute in mixed-type DataFrame
  objects has been implemented (fixes :issue:`135`)
- `read_csv` can read multiple columns into a `MultiIndex`. DataFrame's
  `to_csv` method will properly write out a `MultiIndex` which can be read
  back (:issue:`151`, thanks to Skipper Seabold)
- Wrote fast time series merging / joining methods in Cython. Will be
  integrated later into DataFrame.join and related functions
- Added `ignore_index` option to `DataFrame.append` for combining unindexed
  records stored in a DataFrame

Improvements to existing features
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- Some speed enhancements with internal Index type-checking function
- `DataFrame.rename` has a new `copy` parameter which can rename a DataFrame
  in place
- Enable unstacking by level name (:issue:`142`)
- Enable sortlevel to work by level name (:issue:`141`)
- `read_csv` can automatically "sniff" other kinds of delimiters using
  `csv.Sniffer` (:issue:`146`)
- Improved speed of unit test suite by about 40%
- Exception will not be raised calling `HDFStore.remove` on non-existent node
  with where clause
- Optimized `_ensure_index` function resulting in performance savings in
  type-checking Index objects

API Changes
~~~~~~~~~~~

Bug Fixes
~~~~~~~~~

- Fixed DataFrame constructor bug causing downstream problems (e.g. .copy()
  failing) when passing a Series as the values along with a column name and
  index
- Fixed single-key groupby on DataFrame with as_index=False (:issue:`160`)
- `Series.shift` was failing on integer Series (:issue:`154`)
- `unstack` methods were producing incorrect output in the case of duplicate
  hierarchical labels. An exception will now be raised (:issue:`147`)
- Calling `count` with level argument caused reduceat failure or segfault in
  earlier NumPy (:issue:`169`)
- Fixed `DataFrame.corrwith` to automatically exclude non-numeric data (GH
  :issue:`144`)
- Unicode handling bug fixes in `DataFrame.to_string` (:issue:`138`)
- Excluding OLS degenerate unit test case that was causing platform specific
  failure (:issue:`149`)
- Skip blosc-dependent unit tests for PyTables < 2.2 (:issue:`137`)
- Calling `copy` on `DateRange` did not copy over attributes to the new object
  (:issue:`168`)
- Fix bug in `HDFStore` in which Panel data could be appended to a Table with
  different item order, thus resulting in an incorrect result read back

Thanks
~~~~~~

- Yaroslav Halchenko
- Jeff Reback
- Skipper Seabold
- Dan Lovell
- Nick Pentreath

pandas 0.4.0
------------

**Release date:** 9/12/2011

New Features
~~~~~~~~~~~~

- `pandas.core.sparse` module: "Sparse" (mostly-NA, or some other fill value)
  versions of `Series`, `DataFrame`, and `Panel`. For low-density data, this
  will result in significant performance boosts, and smaller memory
  footprint. Added `to_sparse` methods to `Series`, `DataFrame`, and
  `Panel`. See online documentation for more on these
- Fancy indexing operator on Series / DataFrame, e.g. via .ix operator. Both
  getting and setting of values is supported; however, setting values will only
  currently work on homogeneously-typed DataFrame objects. Things like:

  - series.ix[[d1, d2, d3]]
  - frame.ix[5:10, ['C', 'B', 'A']], frame.ix[5:10, 'A':'C']
  - frame.ix[date1:date2]

- Significantly enhanced `groupby` functionality

  - Can groupby multiple keys, e.g. df.groupby(['key1', 'key2']). Iteration with
    multiple groupings products a flattened tuple
  - "Nuisance" columns (non-aggregatable) will automatically be excluded from
    DataFrame aggregation operations
  - Added automatic "dispatching to Series / DataFrame methods to more easily
    invoke methods on groups. e.g. s.groupby(crit).std() will work even though
    `std` is not implemented on the `GroupBy` class

- Hierarchical / multi-level indexing

  - New the `MultiIndex` class. Integrated `MultiIndex` into `Series` and
    `DataFrame` fancy indexing, slicing, __getitem__ and __setitem,
    reindexing, etc. Added `level` keyword argument to `groupby` to enable
    grouping by a level of a `MultiIndex`

- New data reshaping functions: `stack` and `unstack` on DataFrame and Series

  - Integrate with MultiIndex to enable sophisticated reshaping of data

- `Index` objects (labels for axes) are now capable of holding tuples
- `Series.describe`, `DataFrame.describe`: produces an R-like table of summary
  statistics about each data column
- `DataFrame.quantile`, `Series.quantile` for computing sample quantiles of data
  across requested axis
- Added general `DataFrame.dropna` method to replace `dropIncompleteRows` and
  `dropEmptyRows`, deprecated those.
- `Series` arithmetic methods with optional fill_value for missing data,
  e.g. a.add(b, fill_value=0). If a location is missing for both it will still
  be missing in the result though.
- fill_value option has been added to `DataFrame`.{add, mul, sub, div} methods
  similar to `Series`
- Boolean indexing with `DataFrame` objects: data[data > 0.1] = 0.1 or
  data[data> other] = 1.
- `pytz` / tzinfo support in `DateRange`

  - `tz_localize`, `tz_normalize`, and `tz_validate` methods added

- Added `ExcelFile` class to `pandas.io.parsers` for parsing multiple sheets out
  of a single Excel 2003 document
- `GroupBy` aggregations can now optionally *broadcast*, e.g. produce an object
  of the same size with the aggregated value propagated
- Added `select` function in all data structures: reindex axis based on
  arbitrary criterion (function returning boolean value),
  e.g. frame.select(lambda x: 'foo' in x, axis=1)
- `DataFrame.consolidate` method, API function relating to redesigned internals
- `DataFrame.insert` method for inserting column at a specified location rather
  than the default __setitem__ behavior (which puts it at the end)
- `HDFStore` class in `pandas.io.pytables` has been largely rewritten using
  patches from Jeff Reback from others. It now supports mixed-type `DataFrame`
  and `Series` data and can store `Panel` objects. It also has the option to
  query `DataFrame` and `Panel` data. Loading data from legacy `HDFStore`
  files is supported explicitly in the code
- Added `set_printoptions` method to modify appearance of DataFrame tabular
  output
- `rolling_quantile` functions; a moving version of `Series.quantile` /
  `DataFrame.quantile`
- Generic `rolling_apply` moving window function
- New `drop` method added to `Series`, `DataFrame`, etc. which can drop a set of
  labels from an axis, producing a new object
- `reindex` methods now sport a `copy` option so that data is not forced to be
  copied then the resulting object is indexed the same
- Added `sort_index` methods to Series and Panel. Renamed `DataFrame.sort`
  to `sort_index`. Leaving `DataFrame.sort` for now.
- Added ``skipna`` option to statistical instance methods on all the data
  structures
- `pandas.io.data` module providing a consistent interface for reading time
  series data from several different sources

Improvements to existing features
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- The 2-dimensional `DataFrame` and `DataMatrix` classes have been extensively
  redesigned internally into a single class `DataFrame`, preserving where
  possible their optimal performance characteristics. This should reduce
  confusion from users about which class to use.

  - Note that under the hood there is a new essentially "lazy evaluation"
    scheme within respect to adding columns to DataFrame. During some
    operations, like-typed blocks will be "consolidated" but not before.

- `DataFrame` accessing columns repeatedly is now significantly faster than
  `DataMatrix` used to be in 0.3.0 due to an internal Series caching mechanism
  (which are all views on the underlying data)
- Column ordering for mixed type data is now completely consistent in
  `DataFrame`. In prior releases, there was inconsistent column ordering in
  `DataMatrix`
- Improved console / string formatting of DataMatrix with negative numbers
- Improved tabular data parsing functions, `read_table` and `read_csv`:

  - Added `skiprows` and `na_values` arguments to `pandas.io.parsers` functions
    for more flexible IO
  - `parseCSV` / `read_csv` functions and others in `pandas.io.parsers` now can
    take a list of custom NA values, and also a list of rows to skip

- Can slice `DataFrame` and get a view of the data (when homogeneously typed),
  e.g. frame.xs(idx, copy=False) or frame.ix[idx]
- Many speed optimizations throughout `Series` and `DataFrame`
- Eager evaluation of groups when calling ``groupby`` functions, so if there is
  an exception with the grouping function it will raised immediately versus
  sometime later on when the groups are needed
- `datetools.WeekOfMonth` offset can be parameterized with `n` different than 1
  or -1.
- Statistical methods on DataFrame like `mean`, `std`, `var`, `skew` will now
  ignore non-numerical data. Before a not very useful error message was
  generated. A flag `numeric_only` has been added to `DataFrame.sum` and
  `DataFrame.count` to enable this behavior in those methods if so desired
  (disabled by default)
- `DataFrame.pivot` generalized to enable pivoting multiple columns into a
  `DataFrame` with hierarchical columns
- `DataFrame` constructor can accept structured / record arrays
- `Panel` constructor can accept a dict of DataFrame-like objects. Do not
  need to use `from_dict` anymore (`from_dict` is there to stay, though).

API Changes
~~~~~~~~~~~

- The `DataMatrix` variable now refers to `DataFrame`, will be removed within
  two releases
- `WidePanel` is now known as `Panel`. The `WidePanel` variable in the pandas
  namespace now refers to the renamed `Panel` class
- `LongPanel` and `Panel` / `WidePanel` now no longer have a common
  subclass. `LongPanel` is now a subclass of `DataFrame` having a number of
  additional methods and a hierarchical index instead of the old
  `LongPanelIndex` object, which has been removed. Legacy `LongPanel` pickles
  may not load properly
- Cython is now required to build `pandas` from a development branch. This was
  done to avoid continuing to check in cythonized C files into source
  control. Builds from released source distributions will not require Cython
- Cython code has been moved up to a top level `pandas/src` directory. Cython
  extension modules have been renamed and promoted from the `lib` subpackage to
  the top level, i.e.

  - `pandas.lib.tseries` -> `pandas._tseries`
  - `pandas.lib.sparse` -> `pandas._sparse`

- `DataFrame` pickling format has changed. Backwards compatibility for legacy
  pickles is provided, but it's recommended to consider PyTables-based
  `HDFStore` for storing data with a longer expected shelf life
- A `copy` argument has been added to the `DataFrame` constructor to avoid
  unnecessary copying of data. Data is no longer copied by default when passed
  into the constructor
- Handling of boolean dtype in `DataFrame` has been improved to support storage
  of boolean data with NA / NaN values. Before it was being converted to float64
  so this should not (in theory) cause API breakage
- To optimize performance, Index objects now only check that their labels are
  unique when uniqueness matters (i.e. when someone goes to perform a
  lookup). This is a potentially dangerous tradeoff, but will lead to much
  better performance in many places (like groupby).
- Boolean indexing using Series must now have the same indices (labels)
- Backwards compatibility support for begin/end/nPeriods keyword arguments in
  DateRange class has been removed
- More intuitive / shorter filling aliases `ffill` (for `pad`) and `bfill` (for
  `backfill`) have been added to the functions that use them: `reindex`,
  `asfreq`, `fillna`.
- `pandas.core.mixins` code moved to `pandas.core.generic`
- `buffer` keyword arguments (e.g. `DataFrame.toString`) renamed to `buf` to
  avoid using Python built-in name
- `DataFrame.rows()` removed (use `DataFrame.index`)
- Added deprecation warning to `DataFrame.cols()`, to be removed in next release
- `DataFrame` deprecations and de-camelCasing: `merge`, `asMatrix`,
  `toDataMatrix`, `_firstTimeWithValue`, `_lastTimeWithValue`, `toRecords`,
  `fromRecords`, `tgroupby`, `toString`
- `pandas.io.parsers` method deprecations

  - `parseCSV` is now `read_csv` and keyword arguments have been de-camelCased
  - `parseText` is now `read_table`
  - `parseExcel` is replaced by the `ExcelFile` class and its `parse` method

- `fillMethod` arguments (deprecated in prior release) removed, should be
  replaced with `method`
- `Series.fill`, `DataFrame.fill`, and `Panel.fill` removed, use `fillna`
  instead
- `groupby` functions now exclude NA / NaN values from the list of groups. This
  matches R behavior with NAs in factors e.g. with the `tapply` function
- Removed `parseText`, `parseCSV` and `parseExcel` from pandas namespace
- `Series.combineFunc` renamed to `Series.combine` and made a bit more general
  with a `fill_value` keyword argument defaulting to NaN
- Removed `pandas.core.pytools` module. Code has been moved to
  `pandas.core.common`
- Tacked on `groupName` attribute for groups in GroupBy renamed to `name`
- Panel/LongPanel `dims` attribute renamed to `shape` to be more conforming
- Slicing a `Series` returns a view now
- More Series deprecations / renaming: `toCSV` to `to_csv`, `asOf` to `asof`,
  `merge` to `map`, `applymap` to `apply`, `toDict` to `to_dict`,
  `combineFirst` to `combine_first`. Will print `FutureWarning`.
- `DataFrame.to_csv` does not write an "index" column label by default
  anymore since the output file can be read back without it. However, there
  is a new ``index_label`` argument. So you can do ``index_label='index'`` to
  emulate the old behavior
- `datetools.Week` argument renamed from `dayOfWeek` to `weekday`
- `timeRule` argument in `shift` has been deprecated in favor of using the
  `offset` argument for everything. So you can still pass a time rule string
  to `offset`
- Added optional `encoding` argument to `read_csv`, `read_table`, `to_csv`,
  `from_csv` to handle unicode in Python 2.x

Bug Fixes
~~~~~~~~~

- Column ordering in `pandas.io.parsers.parseCSV` will match CSV in the presence
  of mixed-type data
- Fixed handling of Excel 2003 dates in `pandas.io.parsers`
- `DateRange` caching was happening with high resolution `DateOffset` objects,
  e.g. `DateOffset(seconds=1)`. This has been fixed
- Fixed __truediv__ issue in `DataFrame`
- Fixed `DataFrame.toCSV` bug preventing IO round trips in some cases
- Fixed bug in `Series.plot` causing matplotlib to barf in exceptional cases
- Disabled `Index` objects from being hashable, like ndarrays
- Added `__ne__` implementation to `Index` so that operations like ts[ts != idx]
  will work
- Added `__ne__` implementation to `DataFrame`
- Bug / unintuitive result when calling `fillna` on unordered labels
- Bug calling `sum` on boolean DataFrame
- Bug fix when creating a DataFrame from a dict with scalar values
- Series.{sum, mean, std, ...} now return NA/NaN when the whole Series is NA
- NumPy 1.4 through 1.6 compatibility fixes
- Fixed bug in bias correction in `rolling_cov`, was affecting `rolling_corr`
  too
- R-square value was incorrect in the presence of fixed and time effects in
  the `PanelOLS` classes
- `HDFStore` can handle duplicates in table format, will take

Thanks
~~~~~~

- Joon Ro
- Michael Pennington
- Chris Uga
- Chris Withers
- Jeff Reback
- Ted Square
- Craig Austin
- William Ferreira
- Daniel Fortunov
- Tony Roberts
- Martin Felder
- John Marino
- Tim McNamara
- Justin Berka
- Dieter Vandenbussche
- Shane Conway
- Skipper Seabold
- Chris Jordan-Squire

pandas 0.3.0
------------

**Release date:** February 20, 2011

New features
~~~~~~~~~~~~

- `corrwith` function to compute column- or row-wise correlations between two
  DataFrame objects
- Can boolean-index DataFrame objects, e.g. df[df > 2] = 2, px[px > last_px] = 0
- Added comparison magic methods (__lt__, __gt__, etc.)
- Flexible explicit arithmetic methods (add, mul, sub, div, etc.)
- Added `reindex_like` method
- Added `reindex_like` method to WidePanel
- Convenience functions for accessing SQL-like databases in `pandas.io.sql`
  module
- Added (still experimental) HDFStore class for storing pandas data
  structures using HDF5 / PyTables in `pandas.io.pytables` module
- Added WeekOfMonth date offset
- `pandas.rpy` (experimental) module created, provide some interfacing /
  conversion between rpy2 and pandas

Improvements to existing features
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- Unit test coverage: 100% line coverage of core data structures
- Speed enhancement to rolling_{median, max, min}
- Column ordering between DataFrame and DataMatrix is now consistent: before
  DataFrame would not respect column order
- Improved {Series, DataFrame}.plot methods to be more flexible (can pass
  matplotlib Axis arguments, plot DataFrame columns in multiple subplots,
  etc.)

API Changes
~~~~~~~~~~~

- Exponentially-weighted moment functions in `pandas.stats.moments` have a
  more consistent API and accept a min_periods argument like their regular
  moving counterparts.
- **fillMethod** argument in Series, DataFrame changed to **method**,
  `FutureWarning` added.
- **fill** method in Series, DataFrame/DataMatrix, WidePanel renamed to
  **fillna**, `FutureWarning` added to **fill**
- Renamed **DataFrame.getXS** to **xs**, `FutureWarning` added
- Removed **cap** and **floor** functions from DataFrame, renamed to
  **clip_upper** and **clip_lower** for consistency with NumPy

Bug Fixes
~~~~~~~~~

- Fixed bug in IndexableSkiplist Cython code that was breaking rolling_max
  function
- Numerous numpy.int64-related indexing fixes
- Several NumPy 1.4.0 NaN-handling fixes
- Bug fixes to pandas.io.parsers.parseCSV
- Fixed `DateRange` caching issue with unusual date offsets
- Fixed bug in `DateRange.union`
- Fixed corner case in `IndexableSkiplist` implementation
