#!/usr/bin/env bash


mkdir -p analysis docs data

for my_director in scripts analysis docs data:do
	touch ${my_directory}/README.md
	echo "# ${my_directory}" >> ${my_directory}README.md
done

