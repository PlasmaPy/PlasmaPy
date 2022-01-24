
****************
Code Style Guide
****************


.. important::

   Code is read significantly more often than it is written.


.. tip::

   Every decision we make in writing code involves tradeoffs.



Names
=====

Names are the most fundamental means of communicating the intent and
purpose of code.

* Use named variables instead of named constants.  For example,

* Choose names that are unambiguous, pronounceable, and searchable.

* Avoid purposefully misspelling of names.

* Use only ASCII characters in names that are part of the public API.

* Unicode characters may be used in code that is not user-facing.

.. tip::

   Measure the length of a variable not by the number of characters, but
   instead by the time it takes to understand what it means.

   A variable named ``wx`` would take a long time to understand, so it
   is a long variable name.  A variable named ``nuclear_reaction_energy``
   is short because it can be understood only by reading it.

Functions
=========

Write short functions that do exactly one thing with no side effects.

Units
=====

...
