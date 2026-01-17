## Summary for Opus 4.5

### The Problem
RadioBio survival fraction calculation has two methods (Original O(n²) and Optimized O(n)) that should give identical results but don't.

### Key Observations from User

1. **Offset pattern**: Original and Optimized values look similar if you **shift Original down by one age**
   - Day 34 (pre-injection): Good match with 1-row offset
   - Day 35 (post-injection): Diverge more even with offset
   - See output data in day10.csv

2. **Calling frequency**: `SurvivalProbLookupTableCalc` is called **every hour** (not daily)
   - Takes `currentDay` and `currentHour` as arguments
   - Called inside `Grid.Step()` hourly loop

3. **Age granularity mismatch**:
   - Age is in **days** (integer 1, 2, 3...)
   - But cells born at different **hours** within a day
   - `sliceLength = age * 24` (hours)
   - `birthTime` is in hours

4. **Optimized age 1 always SF=1.0**:
   - New cohorts created with D=0, A=0, G_num=0
   - Need dose history backfilled

5. **Many SF values → 0**: High radiation doses being calculated (frozen tumor so no deaths yet)

### The Code Structure

**HashMap approach (current):**
```java
HashMap<Integer, CohortODEState> cohortStates;  // Key = birthTime in hours
```

**Update pattern:**
1. Every hour: `updateAllCohortStates()` adds current dose to ALL cohorts
2. Then loop calculates SF for ages 1, 2, 3... days
3. Age-to-birthTime: `birthTime = globalTime - (age * 24)`

### Files Available
- `/mnt/project/RadioBio.java` - Current implementation
- User has CSV output showing the offset pattern
- User uploaded images showing the numeric comparison

### The Question
**Should we:**
- Fix the offset issue (seems like age/birthTime indexing problem)?
- Backfill dose history when creating new cohorts?
- Rethink the hourly vs daily update strategy?
- Something else we're missing?

---
