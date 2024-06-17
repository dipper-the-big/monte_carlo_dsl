## 25th june

#### WIP
  1. @02pm: ~~bkl algorithm skeleton.~~

#### TODO
  1. @4pm: ~~skeleton testing with dummy accessory classes.~~

## 26th june

#### WIP
  1. @03pm: ~~Skeleton testing with dummy accessory classes.~~
  2. @03pm: ~~Test for process selection based on dice.~~

#### TODO
  1. @4pm: ~~Skeleton testing with dummy accessory classes.~~

## 18th May

#### Summary
  1. Working okmc system and process collection.

#### TODO
  1. ~~map in System can be changed to unordered_map with a good hash fn. for species.~~
  2. calcCumulative can be made faster if we introduce a speciesCollection that is notified by system of any change in particles and by processes of any change in the rates. I think it will be clear once more things get stable.

#### WIP
  1. ~~Adding Processes for jump and verifying diffusion coefficient.~~

## 22nd Sep

#### TODO
  1. ~~NeighList~~
  2. ~~calcDistMirror~~
  3. ~~caching of recombination radius~~

## 10th Feb

#### TODO
  1. ~~Lazy dirtylist reaction.~~
      * ~~lazy reaction tuples for different types.~~
      * ~~_addDirty with range / iterators.~~
  2. ~~Eager reactions.~~
  3. ~~Fallback reaction.~~

## 16th Feb

#### WIP
  1. Benchmarks and comparison for following problems:
      - Diffusion simulation
      - Kai simulation 
      - ~~Alenya smooth simulation~~
      - He MCD simulation - TODO: Emission wrapper dice cross-check
      - LakiMoca Fe simulation
      - Alenya rough simulation
      - Grain boundary / Dislocation / Structure:

#### TODO
  1. ~~Benchmarks for eager system.~~
  2. ~~Benchmarks for runtime polymorphism with concepts (Sean Parent System).~~
  3. ~~Benchmarks for older systems.~~

## 21st Feb

#### WIP
  1. ~~Correctness check for implemented problems.~~
  2. LakiMoca imlementation.

#### TODO
  1. Implementation for:
     - ~~NeighCellList~~
     - ~~Eager System~~
  2. ~~Benchmarks.~~

## 30th Mar

#### WIP
  1. ~~Alenya example codes for benchmark.~~

#### TODO
  1. ~~Alenya problem benchmark stats.~~
  2. Kai example codes for benchmark.

## 11th April

#### WIP
- Refactored single code of Alenya with:
  1. ~~Basic~~
  2. ~~Cell~~
  3. ~~Eager~~
  4. ~~EagerCell~~
  5. ~~VirtualProcess~~
  6. ~~VirtualProcessEagerCell~~
  7. ~~VirtualReaction~~
  8. ~~VirtualReactionCell~~
  7. ~~VirtualReactionProcess~~
  7. ~~VirtualReactionProcessCell~~
  9. ~~SingleType~~
  10. ~~SingleTypeEagerCell~~
  11. ~~UniquePtr~~
  11. ~~UniquePtrCell~~
  11. ~~UniquePtrEager~~
  12. ~~UniquePtrEagerCell~~
  13. ~~UniquePtrVirtual~~
  14. ~~UniquePtrVirtualEagerCell~~
  15. ~~AllVirtual~~
  16. ~~AllVirtualCell~~
    
#### TODO
- Refactored code for HeMCD, kai and diffCoeff.
- Refactored code for pka aging and one real pka simulation.

## 18th April

#### WIP
- ~~Benchmarking of alenyaAll.~~

## 20th April

#### WIP
- @10AM: ~~Float benchmarks of alenyaAll~~
- @1PM: Refactoring code of HeAll
  1. @2pm : ~~Basic~~
  2. @2pm : ~~Cell~~
  3. @2pm : ~~Eager~~
  4. @3pm : ~~EagerCell~~
  5. @3pm : ~~VirtualProcess~~
  6. @3pm : ~~VirtualProcessEagerCell~~

#### TODO
- @10AM
  - ~~Float benchmarks of alenyaAll~~
  - Benchmark analysis of alenyaAll
  - Refactored code for HeMCD for benchmarks.
  - Reason for eager and neigh interchanging slows down significantly.
  - perf analysis.

## 21st April

#### WIP
- @11AM: Refactoring code of HeAll
  8. @11AM: ~~SingleType~~
  1. @01PM: ~~SingleTypeEagerCell~~
  8. @01PM: ~~VirtualReaction~~
  8. @01PM: ~~VirtualReactionCell~~
  9. @03PM: ~~VirtualReactionProcess~~
  13. @03PM:  ~~UniquePtr~~
  16. @08PM: ~~UniquePtrEagerCell~~

## 23rd April

#### WIP
- @02PM: Refactoring code of HeAll
  17. @02PM: ~~UniquePtrVirtual~~
  18. @04PM: ~~UniquePtrVirtualEagerCell~~
  22. @04PM: ~~AllVirtual~~
  23. @05PM: ~~AllVirtualCell~~
  24. @05PM: ~~NoBufferedRadius~~
- @05PM: ~~Refactoring code of kai problem~~

#### TODO
- @05PM: Refactoring code of HeAll
  25. HybridNeighList
  19. UniquePtrVirtualSingleEagerCell
  20. UniquePtrVirtualSingle
- @11PM: Benchmarking pka problem
- @11PM: Writ

## 13th Sep

#### TODO
- @11am: modernalise and remove compile errors from 'All' examples.