IMAGE=scl3/task_scl_eff_pot_hab


build:
	docker build --no-cache -t $(IMAGE) .

run:
	docker run --rm -it --env-file .env -v `pwd`/src:/app -v `pwd`/.git:/app/.git $(IMAGE) python task.py

shell:
	docker run -it --env-file .env -v `pwd`/src:/app -v `pwd`/.git:/app/.git $(IMAGE) sh

cleanup:
	isort `pwd`/src/*.py
	black `pwd`/src/*.py
	flake8 `pwd`/src/*.py
	mypy `pwd`/src/*.py