#!/bin/bash
python ./../../../src/extra/processResultsInRuntime/watch_and_render.py \
    --config ./template/render_config.toml &
WATCHER_PID=$!

./run.py

# Give the queue time to drain after the last VTK is written
sleep 30

# Then kill the watcher
kill $WATCHER_PID
wait $WATCHER_PID 2>/dev/null
echo "Done."
