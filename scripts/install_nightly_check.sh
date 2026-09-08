#!/usr/bin/env bash
# Install (or remove) the nightly full invariant sweep as a launchd job.
#
#   bash scripts/install_nightly_check.sh            # install, 02:17 nightly
#   bash scripts/install_nightly_check.sh --remove   # uninstall
#
# Runs check_repo_invariants.sh --full, which rebuilds every figure twice and
# once more with data/ hidden. That is minutes of work, which is exactly why it
# belongs on a schedule and not in the pre-push hook.
#
# Log: logs/nightly_check.log (git-ignored), newest run appended with a date
# banner. A failing run leaves a non-zero exit in the log, not a notification --
# check the log, or wire the FAILED line into whatever you already watch.
#
# NOT a Claude agent review. This is the mechanical half. To add the judgment
# half, append a line to the plist command that runs the reviewer headless:
#     claude -p "Use the repo-reviewer agent on the last 24h of commits" \
#       >> logs/nightly_review.log
# That spends tokens unattended, so it is left off by default.
set -uo pipefail
REPO="$(cd "$(dirname "$0")/.." && pwd)"
LABEL="com.stefan.bulkgwa.nightlycheck"
PLIST="$HOME/Library/LaunchAgents/$LABEL.plist"

if [ "${1:-}" = "--remove" ]; then
  launchctl unload "$PLIST" 2>/dev/null || true
  rm -f "$PLIST"
  echo "removed $LABEL"
  exit 0
fi

mkdir -p "$REPO/logs" "$HOME/Library/LaunchAgents"
cat > "$PLIST" <<PLISTEOF
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE plist PUBLIC "-//Apple//DTD PLIST 1.0//EN"
  "http://www.apple.com/DTDs/PropertyList-1.0.dtd">
<plist version="1.0">
<dict>
  <key>Label</key><string>$LABEL</string>
  <key>ProgramArguments</key>
  <array>
    <string>/bin/bash</string>
    <string>-lc</string>
    <string>cd "$REPO" &amp;&amp; { echo; echo "===== \$(date '+%Y-%m-%d %H:%M') ====="; bash scripts/check_repo_invariants.sh --full; echo "exit=\$?"; } >> logs/nightly_check.log 2>&1</string>
  </array>
  <key>StartCalendarInterval</key>
  <dict><key>Hour</key><integer>2</integer><key>Minute</key><integer>17</integer></dict>
  <key>RunAtLoad</key><false/>
</dict>
</plist>
PLISTEOF

launchctl unload "$PLIST" 2>/dev/null || true
launchctl load "$PLIST"
echo "installed $LABEL -- runs 02:17 nightly"
echo "  log:     $REPO/logs/nightly_check.log"
echo "  verify:  launchctl list | grep bulkgwa"
echo "  run now: launchctl start $LABEL"
echo "  remove:  bash scripts/install_nightly_check.sh --remove"
