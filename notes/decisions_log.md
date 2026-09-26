# Decisions log

Judgement calls that do not change the pre-registered plan.

## 2026-09-26

- Git identity for all commits is the owner's global identity, Qasim Hussain, chosen by the owner and set in this repository's local config. No co-author or tool attribution lines are added (rule 0.6).
- The master prompt's opening line says the rules were published on Days 62 to 65, while config.yaml lists days 62 to 66 and rule 0.11 checks days 61 to 66. config.yaml and rule 0.11 are followed.
- All ten download URLs built from the unverified pattern in config.yaml returned HTTP 200 to a header-only request. Script 01 still records each download.
