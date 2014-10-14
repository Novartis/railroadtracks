..
   Copyright 2014 Novartis Institutes for Biomedical Research

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

Persistence
===========

.. note::

   This part of the documentation should only concern users interested in moving the persistence layer
   to a different system to store metadata associated with commands performed.

Tasks to be performed are stored persistently on disk. This is required to ensure that all steps
computed, and the sequence of steps leading to results, are conserved across restarts of the main
process or of the system.

Currently, this is implemented in an SQLite database.


.. automodule:: railroadtracks.hortator
   :members:
