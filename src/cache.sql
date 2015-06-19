-- Copyright 2014-2015 Novartis Institutes for Biomedical Research

-- Licensed under the Apache License, Version 2.0 (the "License");
-- you may not use this file except in compliance with the License.
-- You may obtain a copy of the License at

--     http://www.apache.org/licenses/LICENSE-2.0

-- Unless required by applicable law or agreed to in writing, software
-- distributed under the License is distributed on an "AS IS" BASIS,
-- WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
-- See the License for the specific language governing permissions and
-- limitations under the License.

-- Language-agnostic persistence layer for the DAG of tasks.
 
-- This is implemented as an SQLite database, with the
-- hope that performances will be sufficient. Something fancier
-- will have to be done if this becomes unsufficient.

-- kitchen sink for whatever needs to be stored
-- all values (val) are in JSON
CREATE TABLE misc (
  tag TEXT PRIMARY KEY,
  description TEXT,
  val TEXT -- JSON here
);


-- activities performed by a step, as defined in the model for RNA-Seq
-- (for example: Index, Align, Quantify, Differential Expression)
CREATE TABLE step_activity (
  id INTEGER PRIMARY KEY AUTOINCREMENT,
  label TEXT
);

-- type for a step.
-- a step type can relate to several activities (e.g. Align and Quantify
-- in the same step), or even perform all activities defined in the model
-- (this would be the case of a monolithic pipeline we are bringing in
-- to test)
CREATE TABLE step_type (
  id INTEGER PRIMARY KEY AUTOINCREMENT
);

-- usual design pattern for many-to-many associations
CREATE TABLE step_type2activity (
  step_activity_id INTEGER,
  step_type_id INTEGER,
  FOREIGN KEY(step_activity_id) REFERENCES step_activity(id),
  FOREIGN KEY(step_type_id) REFERENCES step_type(id),
  PRIMARY KEY(step_activity_id, step_type_id)
);

-- variant for a step (that is a class and executable used
-- to define a step - input or output, or parameters are
-- specified through `step_concrete`)
-- Note: executable and version are not factored out into
--       their own table because executable can 
--       toggled to do different activity through the switches
CREATE TABLE step_variant (
  id INTEGER PRIMARY KEY AUTOINCREMENT,
  step_type_id INTEGER,
  executable TEXT,
  cls TEXT, -- class name used to wrap the executable
  version TEXT, -- FIXME: really store it ? it can be inferred from the executable
  FOREIGN KEY(step_type_id) REFERENCES step_type(id)
);

-- Status for a step (for monitoring purposes)
-- Initially the idea is to have 'in-progress' and 'complete'.
-- Should the process running tasks die suddenly,
-- the status 'in-progress' would remain unchanged. A field
-- `last_update` was added to step_concrete to mitigate
-- that (not perfect of completely safe)
CREATE TABLE step_status (
  id INTEGER PRIMARY KEY AUTOINCREMENT,
  label TEXT
);

-- A "concrete" step in a process is a step_variant
-- applied to sources files and yielding target files
-- In addition to that, the parameters and assets are stored in the
-- tables:
-- step_concrete2parameters, step_concrete2srcfile, step_concrete2targetfile
CREATE TABLE step_concrete (
  id INTEGER PRIMARY KEY AUTOINCREMENT,
  step_variant_id INTEGER,
  step_status_id INTEGER,
  time_creation FLOAT, -- system time
  time_t0 FLOAT, -- system time
  time_t1 FLOAT, -- system time
  tag INTEGER NOT NULL,
  FOREIGN KEY(step_variant_id) REFERENCES step_variant(id),
  FOREIGN KEY(step_status_id) REFERENCES step_status(id)
);
CREATE INDEX step_concrete_variant_id_idx ON step_concrete (step_variant_id);

-- Explicitly distinct parameters in the model are labelled.
-- This is used to facilitate the construction of calls to
-- various steps.
-- Labels are determined/defined during the modeling step,
-- and the dependency build graph here does not enforce
-- anything with respect to the model.
-- For example, an aligner will take as input references:
-- - an index file, with the label 'index'
-- - sequencing reads, with the label 'reads'
-- FIXME: since restricting to Illumina, wouldn't 'forward-reads' and 'reverse-reads' be better ?
CREATE TABLE step_parameters (
  id INTEGER PRIMARY KEY AUTOINCREMENT,
  json TEXT UNIQUE
);

-- An 'entity' stored on the disk.
-- This can be a pattern for a file, or group of files.
-- FIXME: because this is based on pattern, a constrain of unicity would be pretty much useless
--        as there can be intersecting patterns.
--        This is a pretty serious design flaw as it is - a move to directories _and_ pattern
--        would be a way to escape this.
-- FIXME: does the extension make any sense if not common to all files (in the case of a group of
--        files ?)
-- FIXME: with this design, each file has _one_ parameter tag.
--        This is working well for the current model (or so I think), and help constrain
--        what can be done and simplify the declarative part and reconstruction of
--        the graph from the database.
--        The big question is whether this constrain would apply to other models.
--        Introducing flexibility here would drive us further toward making a workflow engine <shrug>
CREATE TABLE stored_entity (
  id INTEGER PRIMARY KEY AUTOINCREMENT,
  classname TEXT NOT NULL, -- The Python class name modeling the stored entity
  entityname TEXT NOT NULL, -- this can be NULL (in case associated asset allows None), and the id /might/ be part of the pattern, e.g. ('%(db_id)i.fasta')
  -- step_parameter_id INTEGER,
  -- FOREIGN KEY(step_parameter_id) REFERENCES step_parameter(id)
  UNIQUE(entityname) -- FIXME UNIQUE(classname, entityname)
  -- UNIQUE (pattern_re, extension)
);
CREATE INDEX stored_entity_classname_entityname_idx ON stored_entity (classname, entityname);


-- Special case for sequences
CREATE TABLE stored_sequence (
  id INTEGER PRIMARY KEY,
  classname TEXT NOT NULL -- The Python class name modeling the sequence
);

CREATE TABLE stored_entity2sequence (
  stored_entity_id INTEGER,  -- Stored entity in the sequence
  stored_sequence_id INTEGER,  -- Stored entity in the sequence
  pos INTEGER, -- Position in the sequence
  FOREIGN KEY(stored_entity_id) REFERENCES stored_entity (id),
  FOREIGN KEY(stored_sequence_id) REFERENCES stored_sequence (id),
  UNIQUE(stored_sequence_id, stored_entity_id, pos) -- FIXME primary key ?
);
CREATE INDEX stored_entity2sequence_sequence_idx ON stored_entity2sequence (stored_sequence_id);
CREATE INDEX stored_entity2sequence_entity_idx ON stored_entity2sequence (stored_entity_id);


-- Usual design pattern for many-to-many associations.
-- This permits to navigate the DAG against the direction of its edges
CREATE TABLE step_concrete2srcfile (
  label TEXT,
  step_concrete_id INTEGER NOT NULL,
  stored_entity_id INTEGER,    -- either stored_entity_id or...
  stored_sequence_id INTEGER,  -- ...stored_sequence_id will NOT NULL (see CHECK below)
  --step_parameter_id INTEGER, -- label for the source in the step
  FOREIGN KEY(step_concrete_id) REFERENCES step_concrete(id),
  FOREIGN KEY(stored_entity_id) REFERENCES stored_entity(id),
  FOREIGN KEY(stored_sequence_id) REFERENCES stored_sequence(id),
  --FOREIGN KEY(param_label_id) REFERENCES param_label(id),
  PRIMARY KEY(label, step_concrete_id, stored_entity_id, stored_sequence_id),
  CHECK (((stored_entity_id IS NOT NULL) AND (stored_sequence_id IS NULL))
      OR ((stored_entity_id IS NULL) AND (stored_sequence_id IS NOT NULL))) -- not XOR with SQLite
);

-- Usual design pattern for many-to-many associations.
-- This permits to navigate the DAG in the direction of its edges
CREATE TABLE step_concrete2targetfile (
  label TEXT,
  step_concrete_id INTEGER NOT NULL,
  stored_entity_id INTEGER,    -- either stored_entity_id or...
  stored_sequence_id INTEGER,  -- ...stored_sequence_id will NOT NULL (see CHECK below)
  --step_parameter_id INTEGER, -- name of the target in the step
  FOREIGN KEY(step_concrete_id) REFERENCES step_concrete(id),
  FOREIGN KEY(stored_entity_id) REFERENCES stored_entity(id),
  FOREIGN KEY(stored_sequence_id) REFERENCES stored_sequence(id),
  PRIMARY KEY(step_concrete_id, stored_entity_id, stored_sequence_id),
  CHECK (((stored_entity_id IS NOT NULL) AND (stored_sequence_id IS NULL))
      OR ((stored_entity_id IS NULL) AND (stored_sequence_id IS NOT NULL))) -- not XOR with SQLite
);

CREATE TABLE step_concrete2parameters (
  step_concrete_id INTEGER NOT NULL,
  step_parameters_id INTEGER NOT NULL,
  FOREIGN KEY(step_concrete_id) REFERENCES step_concrete(id),
  FOREIGN KEY(step_parameters_id) REFERENCES step_paramters(id),
  PRIMARY KEY(step_concrete_id, step_parameters_id)
);


